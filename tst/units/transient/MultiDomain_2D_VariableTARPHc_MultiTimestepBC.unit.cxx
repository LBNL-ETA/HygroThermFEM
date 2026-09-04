#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "SlabCreator.hxx"
#include "TestHelpers.hxx"
#include "TestMaterials.hxx"

TEST(MultiDomain_2D_VariableTARPHc_MultiTimestepBC, TestExample_1)
{
    HygroThermFEM::MultiDomain multiDomain;

    const auto & material =
      multiDomain.materials().createSolidMaterial(TestHelper::CottaerSandstone());

    constexpr HygroThermFEM::State state(
      {.temperature = 0.0, .humidity = 0.0, .pressure = 101325.0, .liquidPercent = 1.0});

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
    const std::vector<HygroThermFEM::TARPCoefficients> bcCoeff{{20.0, 0.6},
                                                               {20.0, 0.5},
                                                               {20.0, 0.4},
                                                               {20.0, 0.3},
                                                               {20.0, 0.2},
                                                               {18.0, 0.2},
                                                               {16.0, 0.2},
                                                               {14.0, 0.2},
                                                               {12.0, 0.2},
                                                               {10.0, 0.2}};

    const auto surfaceTilt{90.0};

    multiDomain.createBC_TARPHc(1, 2, bcCoeff, surfaceTilt);

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

    std::vector<std::vector<double>> correctWaterContentSolution{{3.4640534075, 3.4640534075, 0.00336122097966, 0.00336122097967, 6.52286994389e-06, 6.52286994389e-06},
 {5.44169934473, 5.44169934475, 0.0103018353151, 0.0103018353151, 2.97812032171e-05, 2.97812032171e-05},
 {6.48451447967, 6.48451447981, 0.0188394560264, 0.0188394560262, 7.80555109268e-05, 7.80555109269e-05},
 {6.77753301842, 6.77753301868, 0.0282012865604, 0.0282012865602, 0.000157544415391, 0.000157544415391},
 {6.45421147309, 6.45421147347, 0.0377160348728, 0.0377160348724, 0.000271544923395, 0.000271544923396},
 {6.02343378544, 6.02343378562, 0.0471360504168, 0.0471360504165, 0.00042130472766, 0.000421304727661},
 {5.51678214776, 5.51678214805, 0.0562810543505, 0.0562810543501, 0.00060644023091, 0.000606440230912},
 {4.99918367585, 4.99918367611, 0.0647915767321, 0.0647915767317, 0.000824443361255, 0.000824443361257},
 {4.51169878265, 4.51169878299, 0.0724620317653, 0.0724620317652, 0.00107138050555, 0.00107138050555},
 {4.06152005119, 4.06152005152, 0.0793071768127, 0.079307176812, 0.00134285547514, 0.00134285547515}};

    TestHelper::dumpGolden("correctWaterContentSolution", waterContentSolution);
    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{{3.9618296933, 3.9618296933, 2.06012547793, 2.06012547793, 1.54186093011, 1.54186093011},
 {5.53518358448, 5.53518358448, 3.80106543724, 3.80106543724, 3.23272254496, 3.23272254496},
 {6.49387974999, 6.49387974999, 5.13001766381, 5.13001766381, 4.65272886374, 4.65272886374},
 {7.14576829227, 7.14576829226, 6.1198606679, 6.1198606679, 5.75080093255, 5.75080093255},
 {7.58694204143, 7.58694204142, 6.83919288765, 6.83919288765, 6.56542847862, 6.56542847862},
 {7.88322669064, 7.88322669063, 7.35136759833, 7.35136759833, 7.15370842923, 7.15370842923},
 {8.05437573828, 8.05437573827, 7.6963150433, 7.6963150433, 7.55988857377, 7.55988857377},
 {8.12130563069, 8.12130563068, 7.90454710553, 7.90454710553, 7.81793512488, 7.81793512488},
 {8.09609149721, 8.09609149718, 7.99755432664, 7.99755432665, 7.95247628231, 7.95247628231},
 {7.97739099421, 7.97739099419, 7.98552467969, 7.98552467969, 7.97733245, 7.97733245}};

    TestHelper::dumpGolden("correctTemperatureSolution", temperatureSolution);
    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-6);
}
