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

    std::vector<std::vector<double>> correctWaterContentSolution{{5.82829625393, 5.82829625407, 0.00540580634017, 0.00540580634015, 1.0490643726e-05, 1.0490643726e-05},
 {8.55426681031, 8.55426681045, 0.016132736674, 0.0161327366739, 5.0527996886e-05, 5.0527996886e-05},
 {9.2647087277, 9.26470872789, 0.0291841199567, 0.0291841199566, 0.000138040636031, 0.000138040636032},
 {8.47346253134, 8.47346253151, 0.0425078277895, 0.0425078277893, 0.000282709672519, 0.000282709672519},
 {6.50467135305, 6.50467135314, 0.0546411066327, 0.0546411066325, 0.000487524471073, 0.000487524471073},
 {5.03659562085, 5.03659562095, 0.065281544446, 0.065281544446, 0.000748671426611, 0.000748671426612},
 {4.01507830059, 4.0150783005, 0.0737858599486, 0.0737858599469, 0.00105511742358, 0.00105511742358},
 {3.28458692779, 3.28458692791, 0.0807487343467, 0.0807487343448, 0.00139872198926, 0.00139872198927},
 {2.75923556715, 2.75923556703, 0.0865110096024, 0.0865110096012, 0.00177058017718, 0.00177058017719},
 {2.36857058655, 2.36857058646, 0.0913015441361, 0.0913015441348, 0.00216112341915, 0.00216112341916}};

    TestHelper::dumpGolden("correctWaterContentSolution", waterContentSolution);
    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{{6.46596975719, 6.46596975718, 3.36118406544, 3.36118406545, 2.51561191241, 2.51561191241},
 {8.63423773175, 8.63423773174, 5.99429102031, 5.99429102031, 5.11916677108, 5.11916677108},
 {9.66194976388, 9.66194976387, 7.79161364344, 7.79161364344, 7.11933109957, 7.11933109957},
 {10.642773839, 10.6427738389, 9.19229672704, 9.19229672704, 8.6708514202, 8.6708514202},
 {11.1688314364, 11.1688314364, 10.1583745265, 10.1583745265, 9.78423984709, 9.78423984709},
 {11.3536201954, 11.3536201954, 10.7372882785, 10.7372882785, 10.4976428606, 10.4976428606},
 {11.5455478647, 11.5455478647, 11.1315339062, 11.1315339062, 10.9722046866, 10.9722046866},
 {11.5308754677, 11.5308754677, 11.3228980214, 11.3228980214, 11.2348348445, 11.2348348445},
 {11.2893584442, 11.2893584442, 11.2979092901, 11.2979092901, 11.2822203296, 11.2822203296},
 {10.8259239216, 10.8259239216, 11.0539332123, 11.0539332123, 11.111555337, 11.111555337}};

    TestHelper::dumpGolden("correctTemperatureSolution", temperatureSolution);
    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-6);
}
