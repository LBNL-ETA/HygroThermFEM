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

    std::vector<std::vector<double>> correctWaterContentSolution{{8.3440447054, 8.3440447054, 0.00665912097478, 0.00665912097478, 1.29228576309e-05, 1.29228576309e-05},
 {9.57481661031, 9.57481661025, 0.0235678869021, 0.0235678869022, 8.53096752844e-05, 8.53096752843e-05},
 {7.67223907891, 7.67223907887, 0.0430756854524, 0.0430756854525, 0.000272202510658, 0.000272202510658},
 {4.70866152696, 4.70866152694, 0.0599979041521, 0.0599979041522, 0.000610042403431, 0.00061004240343},
 {2.54283103778, 2.54283103802, 0.0714603822579, 0.0714603822605, 0.00112907150947, 0.00112907150947},
 {1.61381523829, 1.61381523837, 0.0805770251789, 0.0805770251821, 0.00188381651966, 0.00188381651965},
 {1.20858375753, 1.2085837576, 0.0873761467716, 0.0873761467746, 0.0028220981627, 0.00282209816267},
 {0.986096666951, 0.986096666997, 0.0928251437848, 0.0928251437873, 0.00390155155245, 0.00390155155241},
 {0.850397712616, 0.850397712641, 0.0972702836631, 0.0972702836654, 0.00507087266875, 0.0050708726687},
 {0.761848577695, 0.761848577706, 0.100917984177, 0.100917984179, 0.00627667187651, 0.00627667187646}};

    TestHelper::dumpGolden("correctWaterContentSolution", waterContentSolution);
    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{{12.2099649291, 12.2099649291, 6.34301191946, 6.34301191946, 4.74729849963, 4.74729849963},
 {16.2922534442, 16.2922534442, 11.3096882462, 11.3096882462, 9.65879514538, 9.65879514538},
 {19.206555786, 19.206555786, 15.2076974201, 15.2076974201, 13.8118135574, 13.8118135574},
 {23.1617646605, 23.1617646606, 19.1696692812, 19.1696692812, 17.821921989, 17.821921989},
 {27.4161316071, 27.416131607, 23.287470473, 23.287470473, 21.9127352203, 21.9127352203},
 {28.1143750371, 28.114375037, 25.623551375, 25.623551375, 24.6904180489, 24.6904180489},
 {28.6562584648, 28.6562584647, 27.0834769809, 27.0834769809, 26.4819895739, 26.4819895739},
 {28.6080484475, 28.6080484474, 27.8020589379, 27.802058938, 27.4706129582, 27.4706129582},
 {28.0236539712, 28.0236539712, 27.8780476283, 27.8780476283, 27.7762669503, 27.7762669503},
 {26.9932922304, 26.9932922304, 27.4084649389, 27.4084649389, 27.5017465034, 27.5017465034}};

    TestHelper::dumpGolden("correctTemperatureSolution", temperatureSolution);
    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-6);
}
