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

    std::vector<std::vector<double>> correctWaterContentSolution{{8.3440432, 8.3440432, 0.00665912024, 0.00665912024, 1.29228562e-05, 1.29228562e-05},
 {9.57270816, 9.57270816, 0.0236334714, 0.0236334714, 8.55217094e-05, 8.55217094e-05},
 {7.66577533, 7.66577533, 0.0431999582, 0.0431999582, 0.000272934272, 0.000272934272},
 {4.70312981, 4.70312981, 0.0601481014, 0.0601481014, 0.000611480873, 0.000611480873},
 {2.53951293, 2.53951293, 0.0716242465, 0.0716242465, 0.00113144662, 0.00113144662},
 {1.6120526, 1.6120526, 0.0807503389, 0.0807503389, 0.00188751167, 0.00188751167},
 {1.20762898, 1.20762898, 0.0875501762, 0.0875501762, 0.00282726125, 0.00282726125},
 {0.985587931, 0.985587931, 0.0929978606, 0.0929978606, 0.00390826239, 0.00390826239},
 {0.850122393, 0.850122393, 0.0974410896, 0.0974410896, 0.00507913808, 0.00507913808},
 {0.761701273, 0.761701273, 0.101086861, 0.101086861, 0.00628644708, 0.00628644708}};

    TestHelper::dumpGolden("correctWaterContentSolution", waterContentSolution);
    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{{12.2156239, 12.2156239, 6.34050025, 6.34050025, 4.74541197, 4.74541197},
 {16.3078312, 16.3078312, 11.3034342, 11.3034342, 9.65359962, 9.65359962},
 {19.2248705, 19.2248705, 15.1981798, 15.1981798, 13.8032728, 13.8032728},
 {23.179048, 23.179048, 19.1590524, 19.1590524, 17.8116146, 17.8116146},
 {27.4285946, 27.4285946, 23.2779364, 23.2779364, 21.9026612, 21.9026612},
 {28.123979, 28.123979, 25.6151266, 25.6151266, 24.6810658, 24.6810658},
 {28.6625959, 28.6625959, 27.0757797, 27.0757797, 26.473227, 26.473227},
 {28.6121123, 28.6121123, 27.7947087, 27.7947087, 27.4621543, 27.4621543},
 {28.0260345, 28.0260345, 27.8709142, 27.8709142, 27.767983, 27.767983},
 {26.9943731, 26.9943731, 27.4015273, 27.4015273, 27.4936311, 27.4936311}};

    TestHelper::dumpGolden("correctTemperatureSolution", temperatureSolution);
    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-6);
}
