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

    std::vector<std::vector<double>> correctWaterContentSolution{{5.82829300489, 5.82829300489, 0.00540581599819, 0.00540581599819, 1.04906624687e-05, 1.04906624687e-05},
 {8.55424768267, 8.55424768267, 0.0161453305933, 0.0161453305933, 5.05631482293e-05, 5.05631482293e-05},
 {9.26467713794, 9.26467713791, 0.0292067498957, 0.0292067498957, 0.000138152352121, 0.000138152352121},
 {8.47342762617, 8.47342762617, 0.0425374113551, 0.0425374113551, 0.000282930350644, 0.000282930350644},
 {6.5046441803, 6.5046441803, 0.0546736323362, 0.0546736323362, 0.000487873904852, 0.000487873904852},
 {5.03658071729, 5.03658071729, 0.0653152303233, 0.0653152303233, 0.000749159874825, 0.000749159874825},
 {4.0150685473, 4.01506854727, 0.0738198794302, 0.0738198794302, 0.00105574891741, 0.00105574891741},
 {3.28458066001, 3.28458066001, 0.0807827719246, 0.0807827719247, 0.00139949881313, 0.00139949881313},
 {2.75923154642, 2.75923154641, 0.0865449346592, 0.0865449346592, 0.0017715024202, 0.0017715024202},
 {2.36856797245, 2.36856797245, 0.09133532404, 0.09133532404, 0.00216218923569, 0.00216218923569}};

    TestHelper::dumpGolden("correctWaterContentSolution", waterContentSolution);
    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{{6.46597006261, 6.46597006261, 3.36118422479, 3.36118422479, 2.51561203167, 2.51561203167},
 {8.63423997272, 8.63423997272, 5.99429218814, 5.99429218814, 5.1191676808, 5.1191676808},
 {9.66195468861, 9.66195468862, 7.79161662722, 7.79161662722, 7.11933358099, 7.11933358099},
 {10.6427815682, 10.6427815682, 9.19230198978, 9.19230198978, 8.67085601851, 8.67085601851},
 {11.1688406448, 11.1688406448, 10.1583816374, 10.1583816374, 9.78424637435, 9.78424637435},
 {11.3536309165, 11.3536309165, 10.7372970836, 10.7372970836, 10.4976511508, 10.4976511508},
 {11.5455581283, 11.5455581283, 11.1315433101, 11.1315433101, 10.9722138741, 10.9722138741},
 {11.530884796, 11.530884796, 11.3229072795, 11.3229072795, 11.2348441523, 11.2348441523},
 {11.289366737, 11.289366737, 11.2979179891, 11.2979179891, 11.2822292515, 11.2822292515},
 {10.8259311694, 10.8259311694, 11.0539411371, 11.0539411371, 11.1115635833, 11.1115635833}};

    TestHelper::dumpGolden("correctTemperatureSolution", temperatureSolution);
    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-6);
}
