#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "SlabCreator.hxx"
#include "TestHelpers.hxx"
#include "TestMaterials.hxx"

TEST(MultiDomain_2D_KimuraHc_MultiTimestepBC, TestExample_1)
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

    // Variable boundary conditions (temperature, wind speed and humidity) over ten timesteps.
    const std::vector<HygroThermFEM::KimuraCoefficients> bcCoeff{
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

    multiDomain.createBC_KimuraHc(1, 2, bcCoeff);

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

    std::vector<std::vector<double>> correctWaterContentSolution{{8.34404319873, 8.34404319873, 0.00665912024366, 0.00665912024366, 1.2922856212e-05, 1.2922856212e-05},
 {9.57390214527, 9.57390214527, 0.0236309884594, 0.0236309884594, 8.55259527985e-05, 8.55259527985e-05},
 {7.66921144249, 7.66921144249, 0.0431925936899, 0.0431925936899, 0.000272994834933, 0.000272994834933},
 {4.70587192503, 4.70587192503, 0.0601415591641, 0.0601415591641, 0.000611735995365, 0.000611735995365},
 {2.540983686, 2.540983686, 0.0716195286015, 0.0716195286015, 0.00113204676087, 0.00113204676087},
 {1.61270472746, 1.61270472746, 0.0807473288888, 0.0807473288888, 0.00188857518216, 0.00188857518216},
 {1.20789462349, 1.20789462349, 0.0875477382003, 0.0875477382003, 0.00282885075392, 0.00282885075392},
 {0.985663179103, 0.985663179104, 0.092995485434, 0.092995485434, 0.00391039275678, 0.00391039275678},
 {0.850113640872, 0.850113640873, 0.0974385083361, 0.0974385083361, 0.00508179603554, 0.00508179603554},
 {0.761656022631, 0.761656022631, 0.101083970561, 0.101083970561, 0.00628959522234, 0.00628959522234}};

    TestHelper::dumpGolden("correctWaterContentSolution", waterContentSolution);
    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{{12.2122939494, 12.2122939494, 6.34281424909, 6.34281424909, 4.74714953559, 4.74714953559},
 {16.2995842899, 16.2995842899, 11.3103653781, 11.3103653781, 9.65925608011, 9.65925608011},
 {19.2161284502, 19.2161284502, 15.2085507492, 15.2085507492, 13.81253948, 13.81253948},
 {23.1720674481, 23.1720674481, 19.1704989319, 19.1704989319, 17.8226601153, 17.8226601153},
 {27.4250615662, 27.4250615662, 23.2885408966, 23.2885408966, 21.9136022997, 21.9136022997},
 {28.1222482449, 28.1222482449, 25.6249919055, 25.6249919055, 24.6915288072, 24.6915288072},
 {28.6625232336, 28.6625232336, 27.0848356874, 27.0848356874, 26.4830427803, 26.4830427803},
 {28.6129754979, 28.6129754979, 27.8031051159, 27.8031051159, 27.471374034, 27.471374034},
 {28.0274487856, 28.0274487856, 27.8787074472, 27.8787074472, 27.7766398427, 27.7766398427},
 {26.9961157068, 26.9961157068, 27.4087382358, 27.4087382358, 27.501726074, 27.501726074}};

    TestHelper::dumpGolden("correctTemperatureSolution", temperatureSolution);
    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-6);
}
