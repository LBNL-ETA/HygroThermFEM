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

    std::vector<std::vector<double>> correctWaterContentSolution{{5.828293, 5.828293, 0.005405816, 0.005405816, 1.04906625e-05, 1.04906625e-05},
 {8.55343377, 8.55343377, 0.0161477857, 0.0161477857, 5.0564131e-05, 5.0564131e-05},
 {9.26202893, 9.26202893, 0.0292149971, 0.0292149971, 0.000138155523, 0.000138155523},
 {8.46860305, 8.46860305, 0.0425512287, 0.0425512287, 0.000282926098, 0.000282926098},
 {6.49901277, 6.49901277, 0.0546896073, 0.0546896073, 0.000487828275, 0.000487828275},
 {5.03224552, 5.03224552, 0.0653279214, 0.0653279214, 0.000749019935, 0.000749019935},
 {4.01184386, 4.01184386, 0.0738298661, 0.0738298661, 0.00105548132, 0.00105548132},
 {3.28233972, 3.28233972, 0.0807904087, 0.0807904087, 0.00139908421, 0.00139908421},
 {2.75774844, 2.75774844, 0.0865508048, 0.0865508048, 0.00177092962, 0.00177092962},
 {2.36763069, 2.36763069, 0.0913399841, 0.0913399841, 0.00216145397, 0.00216145397}};

    TestHelper::dumpGolden("correctWaterContentSolution", waterContentSolution);
    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{{6.47007631, 6.47007631, 3.35942178, 3.35942178, 2.51428791, 2.51428791},
 {8.64324556, 8.64324556, 5.99066601, 5.99066601, 5.11610052, 5.11610052},
 {9.67291715, 9.67291715, 7.786244, 7.786244, 7.11449538, 7.11449538},
 {10.6530035, 10.6530035, 9.18516361, 9.18516361, 8.66421861, 8.66421861},
 {11.1775295, 11.1775295, 10.1502986, 10.1502986, 9.7764148, 9.7764148},
 {11.3613065, 11.3613065, 10.7292199, 10.7292199, 10.4894909, 10.4894909},
 {11.5510891, 11.5510891, 11.1238562, 11.1238562, 10.9642362, 10.9642362},
 {11.5346262, 11.5346262, 11.3156139, 11.3156139, 11.2271856, 11.2271856},
 {11.2916768, 11.2916768, 11.2909875, 11.2909875, 11.2749066, 11.2749066},
 {10.8270891, 10.8270891, 11.0473485, 11.0473485, 11.1045685, 11.1045685}};

    TestHelper::dumpGolden("correctTemperatureSolution", temperatureSolution);
    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-6);
}
