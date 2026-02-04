#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::Nodes;

class MultiDomain_2D_2 : public testing::Test
{
protected:
    void SetUp() override
    {}

    void TearDown() override
    {}
};

TEST_F(MultiDomain_2D_2, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Simple two elements example with moisture and heat transfer.");

    HygroThermFEM::MultiDomain multiDomain;

    // Enter nodes. Arguments are: node number, x-coordinate, y-coordinate

    std::vector<double> gridXCoordinates{0, 0.05, 0.1};

    const double initialTemperature = 0.0;
    const double initialMoistureContent = 0.0;
    const double initialPressure = 101325;
    constexpr auto liquidPercent = 1.0;

    auto state = HygroThermFEM::State(
      initialTemperature, initialMoistureContent, initialPressure, liquidPercent);
    size_t nodeIndex = 0;
    for(auto val : gridXCoordinates)
    {
        ++nodeIndex;
        multiDomain.nodes().createNode(nodeIndex, val, 0.00, state);
        ++nodeIndex;
        multiDomain.nodes().createNode(nodeIndex, val, 0.05, state);
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
        multiDomain.createElement(node1, node2, node3, node4, material.name());
    }

    /// Create Boundary Conditions
    constexpr auto hc = 1.0;
    constexpr auto airTemperature = 20.0;
    constexpr auto humidity = 0.6;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff{airTemperature, hc, humidity};

    multiDomain.createBC_FixedHc(1, 2, bcCoeff);

    constexpr auto dTime = 3600;
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
      {1.089281, 1.089281, 0.001245, 0.001245, 2e-06, 2e-06},
      {2.11223, 2.11223, 0.00378, 0.00378, 1.1e-05, 1.1e-05},
      {3.000406, 3.000406, 0.007521, 0.007521, 3e-05, 3e-05},
      {3.72435, 3.72435, 0.012294, 0.012294, 6.2e-05, 6.2e-05},
      {4.336618, 4.336618, 0.018007, 0.018007, 0.00011, 0.00011},
      {4.866186, 4.866186, 0.024596, 0.024596, 0.00018, 0.00018},
      {5.36041, 5.36041, 0.032008, 0.032008, 0.000273, 0.000273},
      {6.163777, 6.163777, 0.040202, 0.040202, 0.000395, 0.000395},
      {6.88133, 6.88133, 0.049141, 0.049141, 0.000549, 0.000549},
      {7.531928, 7.531928, 0.058814, 0.058814, 0.00074, 0.00074}};

    EXPECT_EQ(waterContentSolution.size(), correctWaterContentSolution.size());

    for(auto i = 0u; i < correctWaterContentSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctWaterContentSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctWaterContentSolution[i][j], waterContentSolution[i][j], 1e-6);
        }
    }

    std::vector<std::vector<double>> correctTemperatureSolution{
      {1.278961, 1.278961, 0.664541, 0.664541, 0.497362, 0.497362},
      {2.103072, 2.103072, 1.389887, 1.389887, 1.165353, 1.165353},
      {2.799476, 2.799476, 2.092455, 2.092455, 1.859221, 1.859221},
      {3.443158, 3.443158, 2.763085, 2.763085, 2.535696, 2.535696},
      {4.052248, 4.052248, 3.402328, 3.402328, 3.184303, 3.184303},
      {4.631945, 4.631945, 4.011701, 4.011701, 3.803543, 3.803543},
      {5.184333, 5.184333, 4.592612, 4.592612, 4.394094, 4.394094},
      {5.710387, 5.710387, 5.146121, 5.146121, 4.956918, 4.956918},
      {6.211405, 6.211405, 5.673432, 5.673432, 5.493160, 5.493160},
      {6.688342, 6.688342, 6.175609, 6.175609, 6.003903, 6.003903}};

    EXPECT_EQ(temperatureSolution.size(), correctTemperatureSolution.size());

    for(auto i = 0u; i < correctTemperatureSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctTemperatureSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctTemperatureSolution[i][j], temperatureSolution[i][j], 1e-6);
        }
    }
}

TEST_F(MultiDomain_2D_2, TestExample_1_Repeat)
{
    SCOPED_TRACE("Begin Test: Repeatability test.");

    HygroThermFEM::MultiDomain multiDomain;

    std::vector<double> gridXCoordinates{0, 0.05, 0.1};

    const double initialTemperature = 0.0;
    const double initialMoistureContent = 0.0;
    const double initialPressure = 101325;
    constexpr auto liquidPercent = 1.0;

    auto state = HygroThermFEM::State(
            initialTemperature, initialMoistureContent, initialPressure, liquidPercent);
    size_t nodeIndex = 0;
    for(auto val : gridXCoordinates)
    {
        ++nodeIndex;
        multiDomain.nodes().createNode(nodeIndex, val, 0.00, state);
        ++nodeIndex;
        multiDomain.nodes().createNode(nodeIndex, val, 0.05, state);
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
        multiDomain.createElement(node1, node2, node3, node4, material.name());
    }

    /// Create Boundary Conditions
    constexpr auto hc = 1.0;
    constexpr auto airTemperature = 20.0;
    constexpr auto humidity = 0.6;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff{airTemperature, hc, humidity};

    multiDomain.createBC_FixedHc(1, 2, bcCoeff);

    constexpr auto dTime = 3600;
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
            {1.089281, 1.089281, 0.001245, 0.001245, 2e-06, 2e-06},
            {2.11223, 2.11223, 0.00378, 0.00378, 1.1e-05, 1.1e-05},
            {3.000406, 3.000406, 0.007521, 0.007521, 3e-05, 3e-05},
            {3.72435, 3.72435, 0.012294, 0.012294, 6.2e-05, 6.2e-05},
            {4.336618, 4.336618, 0.018007, 0.018007, 0.00011, 0.00011},
            {4.866186, 4.866186, 0.024596, 0.024596, 0.00018, 0.00018},
            {5.36041, 5.36041, 0.032008, 0.032008, 0.000273, 0.000273},
            {6.163777, 6.163777, 0.040202, 0.040202, 0.000395, 0.000395},
            {6.88133, 6.88133, 0.049141, 0.049141, 0.000549, 0.000549},
            {7.531928, 7.531928, 0.058814, 0.058814, 0.00074, 0.00074}};

    EXPECT_EQ(waterContentSolution.size(), correctWaterContentSolution.size());

    for(auto i = 0u; i < correctWaterContentSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctWaterContentSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctWaterContentSolution[i][j], waterContentSolution[i][j], 1e-6);
        }
    }

    std::vector<std::vector<double>> correctTemperatureSolution{
            {1.278961, 1.278961, 0.664541, 0.664541, 0.497362, 0.497362},
            {2.103072, 2.103072, 1.389887, 1.389887, 1.165353, 1.165353},
            {2.799476, 2.799476, 2.092455, 2.092455, 1.859221, 1.859221},
            {3.443158, 3.443158, 2.763085, 2.763085, 2.535696, 2.535696},
            {4.052248, 4.052248, 3.402328, 3.402328, 3.184303, 3.184303},
            {4.631945, 4.631945, 4.011701, 4.011701, 3.803543, 3.803543},
            {5.184333, 5.184333, 4.592612, 4.592612, 4.394094, 4.394094},
            {5.710387, 5.710387, 5.146121, 5.146121, 4.956918, 4.956918},
            {6.211405, 6.211405, 5.673432, 5.673432, 5.493160, 5.493160},
            {6.688342, 6.688342, 6.175609, 6.175609, 6.003903, 6.003903}};

    EXPECT_EQ(temperatureSolution.size(), correctTemperatureSolution.size());

    for(auto i = 0u; i < correctTemperatureSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctTemperatureSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctTemperatureSolution[i][j], temperatureSolution[i][j], 1e-6);
        }
    }
}