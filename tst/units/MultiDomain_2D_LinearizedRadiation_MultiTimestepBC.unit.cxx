#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

TEST(MultiDomain_2D_LinearizedRadiation_MultiTimestepBC, TestExample_1)
{
    // Enter nodes. Arguments are: node number, x-coordinate, y-coordinate

    HygroThermFEM::MultiDomain multiDomain;

    std::vector<double> gridXCoordinates{0, 0.05, 0.1};

    const HygroThermFEM::State state({
        .temperature = 0.0,
        .humidity = 0.0,
        .pressure = 101325.0,
        .liquidPercent = 1.0
    });
    for(auto val : gridXCoordinates)
    {
        multiDomain.nodes().createNode({.x = val, .y = 0.00, .state = state});
        multiDomain.nodes().createNode({.x = val, .y = 0.05, .state = state});
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
                        {1, 180}}});

    /// Create elements
    for(size_t i = 1; i <= (multiDomain.nodes().maxIndex() - 2) / 2; ++i)
    {
        const auto node1 = 2u * i + 1u;
        const auto node2 = 2u * i + 2u;
        const auto node3 = 2u * i;
        const auto node4 = 2u * i - 1u;
        multiDomain.createElement({.node1 = node1, .node2 = node2, .node3 = node3, .node4 = node4, .material = material.name()});
    }

    /// Create Boundary Conditions

    // Variable boundary conditions (temperature and humidity) over ten timesteps.
    const std::vector<HygroThermFEM::LinearizedRadiationBCCoefficients> linCoeff{{2.1, 20},
                                                                                 {2.0, 19},
                                                                                 {1.9, 18},
                                                                                 {1.8, 17},
                                                                                 {1.8, 16},
                                                                                 {1.8, 15},
                                                                                 {1.7, 16},
                                                                                 {2.7, 17},
                                                                                 {3.7, 18},
                                                                                 {4.7, 19}};

    multiDomain.createBC_LinearizedRadiation(1, 2, linCoeff);

    constexpr auto dTime = 3600;
    constexpr auto nSteps = 10;

    auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
    auto humidities = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);
    std::vector<std::vector<double>> temperatureSolution;
    std::vector<std::vector<double>> waterContentSolution;
    size_t timestepIndex{0u};

    for(auto i = 0; i < nSteps; ++i)
    {
        auto aSolution = multiDomain.transient(temperatures, humidities, dTime, timestepIndex);
        temperatureSolution.push_back(aSolution.temperature);
        waterContentSolution.push_back(aSolution.waterContent);
        temperatures = aSolution.temperature;
        humidities = aSolution.humidity;
        ++timestepIndex;
    }

    std::vector<std::vector<double>> correctWaterContentSolution{{0, 0, 0, 0, 0, 0},
                                                                 {0, 0, 0, 0, 0, 0},
                                                                 {0, 0, 0, 0, 0, 0},
                                                                 {0, 0, 0, 0, 0, 0},
                                                                 {0, 0, 0, 0, 0, 0},
                                                                 {0, 0, 0, 0, 0, 0},
                                                                 {0, 0, 0, 0, 0, 0},
                                                                 {0, 0, 0, 0, 0, 0},
                                                                 {0, 0, 0, 0, 0, 0},
                                                                 {0, 0, 0, 0, 0, 0}};

    EXPECT_EQ(waterContentSolution.size(), correctWaterContentSolution.size());

    for(auto i = 0u; i < correctWaterContentSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctWaterContentSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctWaterContentSolution[i][j], waterContentSolution[i][j], 1e-6);
        }
    }

    std::vector<std::vector<double>> correctTemperatureSolution{
      {1.333840, 1.333840, 0.693325, 0.693325, 0.518905, 0.518905},
      {2.080571, 2.080571, 1.391602, 1.391602, 1.172058, 1.172058},
      {2.633593, 2.633593, 2.008475, 2.008475, 1.798058, 1.798058},
      {3.084483, 3.084483, 2.540265, 2.540265, 2.353548, 2.353548},
      {3.500055, 3.500055, 3.014744, 3.014744, 2.848407, 2.848407},
      {3.869656, 3.869656, 3.437373, 3.437373, 3.289206, 3.289206},
      {4.270165, 4.270165, 3.850880, 3.850880, 3.709580, 3.709580},
      {5.094680, 5.094680, 4.478925, 4.478925, 4.285381, 4.285381},
      {6.188544, 6.188544, 5.342270, 5.342270, 5.076388, 5.076388},
      {7.489070, 7.489070, 6.423400, 6.423400, 6.084532, 6.084532}};

    EXPECT_EQ(temperatureSolution.size(), correctTemperatureSolution.size());

    for(auto i = 0u; i < correctTemperatureSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctTemperatureSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctTemperatureSolution[i][j], temperatureSolution[i][j], 1e-6);
        }
    }
}
