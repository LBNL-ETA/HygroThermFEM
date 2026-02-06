#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "SlabCreator.hxx"
#include "TestMaterials.hxx"

TEST(MultiDomain_2D_BlackBody_MultiTimestepBC, TestExample_1)
{
    HygroThermFEM::MultiDomain multiDomain;

    const auto & material = multiDomain.materials().createSolidMaterial(TestHelper::CottaerSandstone());

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
    const std::vector<HygroThermFEM::BlackBodyRadiationBCCoefficients> bcBlackBody{{0.8, 20},
                                                                                   {0.9, 19},
                                                                                   {0.9, 18},
                                                                                   {0.9, 17},
                                                                                   {0.9, 16},
                                                                                   {0.8, 15},
                                                                                   {0.7, 16},
                                                                                   {0.7, 17},
                                                                                   {0.7, 18},
                                                                                   {0.7, 19}};

    multiDomain.createBC_BlackBodyRadiation(1, 2, bcBlackBody);

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
      {2.489003, 2.489003, 1.293774, 1.293774, 0.968299, 0.968299},
      {4.120427, 4.120427, 2.720497, 2.720497, 2.279696, 2.279696},
      {5.269488, 5.269488, 3.987812, 3.987812, 3.558101, 3.558101},
      {6.181914, 6.181914, 5.072106, 5.072106, 4.691228, 4.691228},
      {6.922008, 6.922008, 5.983873, 5.983873, 5.658682, 5.658682},
      {7.397349, 7.397349, 6.676069, 6.676069, 6.420125, 6.420125},
      {7.941282, 7.941282, 7.300253, 7.300253, 7.078839, 7.078839},
      {8.591005, 8.591005, 7.942228, 7.942228, 7.725025, 7.725025},
      {9.285293, 9.285293, 8.611945, 8.611945, 8.388823, 8.388823},
      {10.009687, 10.009687, 9.30931, 9.30931, 9.077743, 9.077743}};

    EXPECT_EQ(temperatureSolution.size(), correctTemperatureSolution.size());

    for(auto i = 0u; i < correctTemperatureSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctTemperatureSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctTemperatureSolution[i][j], temperatureSolution[i][j], 1e-6);
        }
    }
}
