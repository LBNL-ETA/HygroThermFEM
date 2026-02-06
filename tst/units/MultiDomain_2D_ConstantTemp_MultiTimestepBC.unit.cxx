#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "SlabCreator.hxx"
#include "TestMaterials.hxx"

TEST(MultiDomain_2D_ConstantTemp_MultiTimestepBC, TestExample_1)
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
    const std::vector<double> bcTemperatures{20, 19, 18, 17, 16, 15, 16, 17, 18, 19};

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
