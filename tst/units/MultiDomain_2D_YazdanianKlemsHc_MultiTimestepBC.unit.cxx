#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "SlabCreator.hxx"
#include "TestHelpers.hxx"
#include "TestMaterials.hxx"

TEST(MultiDomain_2D_YazdanianKlemsHc_MultiTimestepBC, TestExample_1)
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

    using HygroThermFEM::WindDirection;

    // Variable boundary conditions (temperature, humidity, wind speed and wind direction) over ten timesteps.
    const std::vector<HygroThermFEM::YazdanianKlemsCoefficients> bcCoeff{
      {20.0, 0.6, 3.0, WindDirection::Windward},
      {20.0, 0.5, 3.0, WindDirection::Windward},
      {20.0, 0.4, 3.0, WindDirection::Windward},
      {20.0, 0.3, 4.0, WindDirection::Windward},
      {20.0, 0.2, 4.2, WindDirection::Windward},
      {18.0, 0.2, 4.6, WindDirection::Leeward},
      {16.0, 0.2, 5.0, WindDirection::Leeward},
      {14.0, 0.2, 5.3, WindDirection::Leeward},
      {12.0, 0.2, 5.5, WindDirection::Leeward},
      {10.0, 0.2, 5.9, WindDirection::Leeward}};

    multiDomain.createBC_YazdanianKlemsHc(1, 2, bcCoeff);

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
      {4.421910, 4.421910, 0.004642, 0.004642, 0.000010, 0.000010},
      {6.948913, 6.948913, 0.012497, 0.012497, 0.000043, 0.000043},
      {7.969250, 7.969250, 0.022601, 0.022601, 0.000113, 0.000113},
      {7.453109, 7.453109, 0.033553, 0.033553, 0.000232, 0.000232},
      {5.509753, 5.509753, 0.043839, 0.043839, 0.000405, 0.000405},
      {4.576782, 4.576782, 0.053225, 0.053225, 0.000630, 0.000630},
      {3.827251, 3.827251, 0.061444, 0.061444, 0.000901, 0.000901},
      {3.192876, 3.192876, 0.068485, 0.068485, 0.001212, 0.001212},
      {2.689526, 2.689526, 0.074463, 0.074463, 0.001554, 0.001554},
      {2.313827, 2.313827, 0.079519, 0.079519, 0.001917, 0.001917}};

    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{
      {6.671904, 6.671904, 3.465066, 3.465066, 2.593355, 2.593355},
      {8.943530, 8.943530, 6.194599, 6.194599, 5.288623, 5.288623},
      {10.002743, 10.002743, 8.051797, 8.051797, 7.356649, 7.356649},
      {11.017891, 11.017891, 9.499455, 9.499455, 8.960367, 8.960367},
      {11.588020, 11.588020, 10.512348, 10.512348, 10.121890, 10.121890},
      {11.787801, 11.787801, 11.122641, 11.122641, 10.870855, 10.870855},
      {11.902599, 11.902599, 11.493919, 11.493919, 11.337146, 11.337146},
      {11.826965, 11.826965, 11.645643, 11.645643, 11.568006, 11.568006},
      {11.546816, 11.546816, 11.583489, 11.583489, 11.579567, 11.579567},
      {11.049692, 11.049692, 11.305091, 11.305091, 11.374116, 11.374116}};

    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-6);
}
