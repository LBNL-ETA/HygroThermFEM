#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::Nodes;
using HygroThermFEM::State;

////////////////////////////////////////////////////////////////////////////////////////////////////
/// This is simple two elements multi-domain example without boundary conditions. Initial
/// temperature and moisture distribution is not same in every node. This case should prove
/// that domain will try to reach equilibrium
////////////////////////////////////////////////////////////////////////////////////////////////////

TEST(MultiDomain_2D_1, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Simple two elements example with moisture and heat transfer.");

    HygroThermFEM::MultiDomain multiDomain;

    // Enter nodes. Arguments are: node number, x-coordinate, y-coordinate

    std::vector<double> gridXCoordinates{0, 0.05, 0.1};

    const double initialTemperature = 20.0;
    const double initialMoistureContent = 0.0;
    const double initialPressure = 101325.0;

    auto tVal = 0.0;
    auto deltaT = 10.0;
    auto hVal = 0.0;
    auto deltaH = 0.1;
    for(auto val : gridXCoordinates)
    {
        multiDomain.nodes().createNode(
          {.x = val,
           .y = 0.00,
           .state = State{
              .temperature = initialTemperature + tVal,
              .humidity = initialMoistureContent + hVal,
              .pressure = initialPressure,
              .liquidPercent = 0
          }});
        multiDomain.nodes().createNode(
          {.x = val,
           .y = 0.05,
           .state = State{
              .temperature = initialTemperature + tVal,
              .humidity = initialMoistureContent + hVal,
              .pressure = initialPressure,
              .liquidPercent = 0
          }});
        tVal += deltaT;
        hVal += deltaH;
    }

    // Material Properties (Cottaer Sandstone)
    const auto & material = multiDomain.materials().createSolidMaterial({
        .name = "Cottaer Sandstone",
        .thermalConductivityDry = 1.8,
        .density = 2050.0,
        .porosity = 0.22,
        .heatCapacity = 850.0,
        .diffusionResistanceFactor = 15.0,
        .thermalConductivityMoistureDependent = {{0.0, 1.8}, {180, 1.8}},
        .moistureDependentMeasurementTemperature = 0,
        .thermalConductivityTemperatureDependent = {{0.0, 1.8}, {1, 1.8}},
        .temperatureDependentMeasurementHumidity = 0,
        .liquidTransportCurve = {{0, 0},
                                 {27, 1E-8},
                                 {45, 1.1E-8},
                                 {90, 2E-8},
                                 {126, 3.5E-8},
                                 {144, 5E-8},
                                 {162, 1E-7},
                                 {171, 2E-7},
                                 {180, 7E-7}},
        .sorptionCurve = {{0, 0},
                          {0.5, 5.3},
                          {0.65, 8.4},
                          {0.8, 12},
                          {0.93, 17},
                          {0.95, 25},
                          {0.99, 63},
                          {0.995, 83},
                          {0.999, 120},
                          {1, 180}}
    });

    /// Create elements
    for(size_t i = 1; i <= (multiDomain.nodes().maxIndex() - 2) / 2; ++i)
    {
        const auto node1 = 2u * i + 1u;
        const auto node2 = 2u * i + 2u;
        const auto node3 = 2u * i;
        const auto node4 = 2u * i - 1u;
        multiDomain.createElement({.node1 = node1, .node2 = node2, .node3 = node3, .node4 = node4, .material = material.name()});
    }

    constexpr auto dTime = 360;
    constexpr auto nSteps = 10;

    auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
    auto humidities = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);
    std::vector<std::vector<double>> temperatureSolution;
    std::vector<std::vector<double>> waterContentSolution;

    for(auto i = 0; i < nSteps; ++i)
    {
        auto aSolution = multiDomain.transient(temperatures, humidities, dTime);
        temperatureSolution.push_back(aSolution.temperature);
        waterContentSolution.push_back(aSolution.waterContent);
        temperatures = aSolution.temperature;
        humidities = aSolution.humidity;
    }

    std::vector<std::vector<double>> correctWaterContentSolution{
      {0.000944, 0.000944, 1.060053, 1.060053, 2.117104, 2.117104},
      {0.001987, 0.001987, 1.060083, 1.060083, 2.114588, 2.114588},
      {0.003113, 0.003113, 1.060099, 1.060099, 2.112348, 2.112348},
      {0.004305, 0.004305, 1.060107, 1.060107, 2.110308, 2.110308},
      {0.005551, 0.005551, 1.060109, 1.060109, 2.108417, 2.108417},
      {0.006839, 0.006839, 1.060109, 1.060109, 2.106638, 2.106638},
      {0.008159, 0.008159, 1.060107, 1.060107, 2.104943, 2.104943},
      {0.009505, 0.009505, 1.060103, 1.060103, 2.103312, 2.103312},
      {0.01087, 0.01087, 1.060099, 1.060099, 2.101731, 2.101731},
      {0.01225, 0.01225, 1.060094, 1.060094, 2.100187, 2.100187}};

    EXPECT_EQ(waterContentSolution.size(), correctWaterContentSolution.size());

    for(auto i = 0u; i < correctWaterContentSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctWaterContentSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctWaterContentSolution[i][j], waterContentSolution[i][j], 1e-6);
        }
    }

    std::vector<std::vector<double>> correctTemperatureSolution{
      {22.291791, 22.291791, 30.001060, 30.001060, 37.718922, 37.718922},
      {24.058697, 24.058697, 30.002681, 30.002681, 35.958658, 35.958658},
      {25.421047, 25.421047, 30.004433, 30.004433, 34.600417, 34.600417},
      {26.471542, 26.471542, 30.006097, 30.006097, 33.552457, 33.552457},
      {27.281613, 27.281613, 30.007575, 30.007575, 32.743944, 32.743944},
      {27.906314, 27.906314, 30.008837, 30.008837, 32.120196, 32.120196},
      {28.388083, 28.388083, 30.009887, 30.009887, 31.639007, 31.639007},
      {28.759632, 28.759632, 30.010743, 30.010743, 31.267808, 31.267808},
      {29.046186, 29.046186, 30.011433, 30.011433, 30.981463, 30.981463},
      {29.267192, 29.267192, 30.011984, 30.011984, 30.760580, 30.760580}};

    EXPECT_EQ(temperatureSolution.size(), correctTemperatureSolution.size());

    for(auto i = 0u; i < correctTemperatureSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctTemperatureSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctTemperatureSolution[i][j], temperatureSolution[i][j], 1e-6);
        }
    }
}

TEST(MultiDomain_2D_1, TestExample_1_Repeat)
{
    SCOPED_TRACE("Begin Test: Repeatability test.");

    HygroThermFEM::MultiDomain multiDomain;

    // Enter nodes. Arguments are: node number, x-coordinate, y-coordinate

    std::vector<double> gridXCoordinates{0, 0.05, 0.1};

    const double initialTemperature = 20.0;
    const double initialMoistureContent = 0.0;
    const double initialPressure = 101325.0;

    auto tVal = 0.0;
    auto deltaT = 10.0;
    auto hVal = 0.0;
    auto deltaH = 0.1;
    for(auto val : gridXCoordinates)
    {
        multiDomain.nodes().createNode(
                {.x = val,
                 .y = 0.00,
                 .state = State{
                    .temperature = initialTemperature + tVal,
                    .humidity = initialMoistureContent + hVal,
                    .pressure = initialPressure,
                    .liquidPercent = 0
                }});
        multiDomain.nodes().createNode(
                {.x = val,
                 .y = 0.05,
                 .state = State{
                    .temperature = initialTemperature + tVal,
                    .humidity = initialMoistureContent + hVal,
                    .pressure = initialPressure,
                    .liquidPercent = 0
                }});
        tVal += deltaT;
        hVal += deltaH;
    }

    // Material Properties (Cottaer Sandstone)
    const auto & material = multiDomain.materials().createSolidMaterial({
        .name = "Cottaer Sandstone",
        .thermalConductivityDry = 1.8,
        .density = 2050.0,
        .porosity = 0.22,
        .heatCapacity = 850.0,
        .diffusionResistanceFactor = 15.0,
        .thermalConductivityMoistureDependent = {{0.0, 1.8}, {180, 1.8}},
        .moistureDependentMeasurementTemperature = 0,
        .thermalConductivityTemperatureDependent = {{0.0, 1.8}, {1, 1.8}},
        .temperatureDependentMeasurementHumidity = 0,
        .liquidTransportCurve = {{0, 0},
                                 {27, 1E-8},
                                 {45, 1.1E-8},
                                 {90, 2E-8},
                                 {126, 3.5E-8},
                                 {144, 5E-8},
                                 {162, 1E-7},
                                 {171, 2E-7},
                                 {180, 7E-7}},
        .sorptionCurve = {{0, 0},
                          {0.5, 5.3},
                          {0.65, 8.4},
                          {0.8, 12},
                          {0.93, 17},
                          {0.95, 25},
                          {0.99, 63},
                          {0.995, 83},
                          {0.999, 120},
                          {1, 180}}
    });

    /// Create elements
    for(size_t i = 1; i <= (multiDomain.nodes().maxIndex() - 2) / 2; ++i)
    {
        const auto node1 = 2u * i + 1u;
        const auto node2 = 2u * i + 2u;
        const auto node3 = 2u * i;
        const auto node4 = 2u * i - 1u;
        multiDomain.createElement({.node1 = node1, .node2 = node2, .node3 = node3, .node4 = node4, .material = material.name()});
    }

    constexpr auto dTime = 360;
    constexpr auto nSteps = 10;

    auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
    auto humidities = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);
    std::vector<std::vector<double>> temperatureSolution;
    std::vector<std::vector<double>> waterContentSolution;

    for(auto i = 0; i < nSteps; ++i)
    {
        auto aSolution = multiDomain.transient(temperatures, humidities, dTime);
        temperatureSolution.push_back(aSolution.temperature);
        waterContentSolution.push_back(aSolution.waterContent);
        temperatures = aSolution.temperature;
        humidities = aSolution.humidity;
    }

    std::vector<std::vector<double>> correctWaterContentSolution{
            {0.000944, 0.000944, 1.060053, 1.060053, 2.117104, 2.117104},
            {0.001987, 0.001987, 1.060083, 1.060083, 2.114588, 2.114588},
            {0.003113, 0.003113, 1.060099, 1.060099, 2.112348, 2.112348},
            {0.004305, 0.004305, 1.060107, 1.060107, 2.110308, 2.110308},
            {0.005551, 0.005551, 1.060109, 1.060109, 2.108417, 2.108417},
            {0.006839, 0.006839, 1.060109, 1.060109, 2.106638, 2.106638},
            {0.008159, 0.008159, 1.060107, 1.060107, 2.104943, 2.104943},
            {0.009505, 0.009505, 1.060103, 1.060103, 2.103312, 2.103312},
            {0.01087, 0.01087, 1.060099, 1.060099, 2.101731, 2.101731},
            {0.01225, 0.01225, 1.060094, 1.060094, 2.100187, 2.100187}};

    EXPECT_EQ(waterContentSolution.size(), correctWaterContentSolution.size());

    for(auto i = 0u; i < correctWaterContentSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctWaterContentSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctWaterContentSolution[i][j], waterContentSolution[i][j], 1e-6);
        }
    }

    std::vector<std::vector<double>> correctTemperatureSolution{
            {22.291791, 22.291791, 30.001060, 30.001060, 37.718922, 37.718922},
            {24.058697, 24.058697, 30.002681, 30.002681, 35.958658, 35.958658},
            {25.421047, 25.421047, 30.004433, 30.004433, 34.600417, 34.600417},
            {26.471542, 26.471542, 30.006097, 30.006097, 33.552457, 33.552457},
            {27.281613, 27.281613, 30.007575, 30.007575, 32.743944, 32.743944},
            {27.906314, 27.906314, 30.008837, 30.008837, 32.120196, 32.120196},
            {28.388083, 28.388083, 30.009887, 30.009887, 31.639007, 31.639007},
            {28.759632, 28.759632, 30.010743, 30.010743, 31.267808, 31.267808},
            {29.046186, 29.046186, 30.011433, 30.011433, 30.981463, 30.981463},
            {29.267192, 29.267192, 30.011984, 30.011984, 30.760580, 30.760580}};

    EXPECT_EQ(temperatureSolution.size(), correctTemperatureSolution.size());

    for(auto i = 0u; i < correctTemperatureSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctTemperatureSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctTemperatureSolution[i][j], temperatureSolution[i][j], 1e-6);
        }
    }
}