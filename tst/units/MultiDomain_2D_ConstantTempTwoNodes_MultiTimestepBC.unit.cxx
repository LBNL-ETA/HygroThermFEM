#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "SlabCreator.hxx"

TEST(MultiDomain_2D_ConstantTempTwoNodes_MultiTimestepBC, TestExample_1)
{
    HygroThermFEM::MultiDomain multiDomain;

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

    constexpr HygroThermFEM::State state({
        .temperature = 0.0,
        .humidity = 0.0,
        .pressure = 101325.0,
        .liquidPercent = 1.0
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

    // Variable boundary conditions (temperature and humidity) over ten timesteps.
    const std::vector<HygroThermFEM::ConstantBCTemperatures> bcTemperatures{{20, 20},
                                                                            {19, 19},
                                                                            {18, 18},
                                                                            {17, 17},
                                                                            {16, 16},
                                                                            {15, 15},
                                                                            {16, 16},
                                                                            {17, 17},
                                                                            {18, 18},
                                                                            {19, 19}};

    multiDomain.createBC_FixedTemperature(1, 2, bcTemperatures);

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
      {20, 20, 10.395924, 10.395924, 7.7806210, 7.7806210},
      {19, 19, 14.526299, 14.526299, 12.829288, 12.829288},
      {18, 18, 16.110006, 16.110006, 15.284675, 15.284675},
      {17, 17, 16.464697, 16.464697, 16.167839, 16.167839},
      {16, 16, 16.184330, 16.184330, 16.180181, 16.180181},
      {15, 15, 15.568177, 15.568177, 15.722139, 15.722139},
      {16, 16, 15.812770, 15.812770, 15.789970, 15.789970},
      {17, 17, 16.426906, 16.426906, 16.266672, 16.266672},
      {18, 18, 17.223641, 17.223641, 16.982896, 16.982896},
      {19, 19, 18.115504, 18.115504, 17.830574, 17.830574}};

    EXPECT_EQ(temperatureSolution.size(), correctTemperatureSolution.size());

    for(auto i = 0u; i < correctTemperatureSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctTemperatureSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctTemperatureSolution[i][j], temperatureSolution[i][j], 1e-6);
        }
    }
}
