#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "SlabCreator.hxx"

TEST(MultiDomain_2D_ASHRAEOutsideHc_MultiTimestepBC, TestExample_1)
{
    HygroThermFEM::MultiDomain multiDomain;

    const HygroThermFEM::State state({
        .temperature = 0.0,
        .humidity = 0.0,
        .pressure = 101325.0,
        .liquidPercent = 1.0
    });

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

    TestHelper::SlabBuilder(multiDomain)
        .gridXCoordinates({0, 0.05, 0.1})
        .height(0.05)
        .material(material.name())
        .state(state)
        .startCorner(TestHelper::StartCorner::BottomRight)
        .direction(TestHelper::Direction::CounterClockwise)
        .build();

    /// Create Boundary Conditions

    // Variable boundary conditions (temperature, humidity and wind speed) over ten timesteps.
    const std::vector<HygroThermFEM::ASHRAEOutsideCoefficients> bcCoeff{{20.0, 0.6, 3},
                                                                       {20.0, 0.5, 3},
                                                                       {20.0, 0.4, 3},
                                                                       {20.0, 0.3, 4},
                                                                       {20.0, 0.2, 4.2},
                                                                       {18.0, 0.2, 4.6},
                                                                       {16.0, 0.2, 5},
                                                                       {14.0, 0.2, 5.3},
                                                                       {12.0, 0.2, 5.5},
                                                                       {10.0, 0.2, 5.9}};

    multiDomain.createBC_ASHRAEOutsideHc(1, 2, bcCoeff);

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

    std::vector<std::vector<double>> correctWaterContentSolution{
      {5.839472,5.839472,0.008244,0.008244,2.3e-05,2.3e-05},
      {8.360858,8.360858,0.020737,0.020737,9.9e-05,9.9e-05},
      {7.723784,7.723784,0.03439,0.03439,0.000248,0.000248},
      {5.210132,5.210132,0.046633,0.046633,0.000473,0.000473},
      {3.487302,3.487302,0.055718,0.055718,0.000765,0.000765},
      {2.673549,2.673549,0.06306,0.06306,0.001115,0.001115},
      {2.299216,2.299216,0.069381,0.069381,0.001509,0.001509},
      {2.098217,2.098217,0.074977,0.074977,0.001929,0.001929},
      {1.975002,1.975002,0.079962,0.079962,0.00236,0.00236},
      {1.895037,1.895037,0.084394,0.084394,0.002787,0.002787}};

    EXPECT_EQ(waterContentSolution.size(), correctWaterContentSolution.size());

    for(auto i = 0u; i < correctWaterContentSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctWaterContentSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctWaterContentSolution[i][j], waterContentSolution[i][j], 1e-6);
        }
    }

    std::vector<std::vector<double>> correctTemperatureSolution{
      {10.519153,10.519153,5.461413,5.461413,4.087475,4.087475},
      {12.721213,12.721213,9.048704,9.048704,7.800583,7.800583},
      {13.533839,13.533839,11.212284,11.212284,10.353969,10.353969},
      {14.527589,14.527589,12.820343,12.820343,12.199838,12.199838},
      {15.451702,15.451702,14.105056,14.105056,13.625717,13.625717},
      {15.79481,15.79481,14.919424,14.919424,14.593923,14.593923},
      {15.460193,15.460193,15.157125,15.157125,15.015402,15.015402},
      {14.620575,14.620575,14.859221,14.859221,14.898481,14.898481},
      {13.435955,13.435955,14.124344,14.124344,14.319072,14.319072},
      {11.978455,11.978455,13.034402,13.034402,13.357576,13.357576}};

    EXPECT_EQ(temperatureSolution.size(), correctTemperatureSolution.size());

    for(auto i = 0u; i < correctTemperatureSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctTemperatureSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctTemperatureSolution[i][j], temperatureSolution[i][j], 1e-6);
        }
    }
}
