#include <gtest/gtest.h>
#include <set>

#include "HygroThermFEM2D.hxx"
#include "TestMaterials.hxx"
#include "TestHelpers.hxx"

using HygroThermFEM::State;

TEST(TestModelWithFrameCavity3, TestSingleFrameCavity)
{
    SCOPED_TRACE("Begin Test: Model with single frame cavity.");

    // Create multi-domain before materials
    HygroThermFEM::MultiDomain multiDomain;

    std::vector<double> gridX{0, 0.01, 0.02, 0.03, 0.04};
    std::vector<double> gridY{0, 0.05, 0.1, 0.15, 0.2};

    constexpr auto initialTemperature{20.0};
    constexpr auto initialHumidity{0.0};
    constexpr auto initialPressure{101325.0};

    constexpr State state({
        .temperature = initialTemperature,
        .humidity = initialHumidity,
        .pressure = initialPressure
    });

    // Crating grid nodes
    for(auto yVal : gridY)
    {
        for(auto xVal : gridX)
        {
            multiDomain.nodes().createNode({.x = xVal, .y = yVal, .state = state});
        }
    }

    auto params = TestHelper::CottaerSandstone();
    params.name = "Material 1";
    params.liquidTransportCurve = {{0, 0}, {180, 2e-6}};
    params.sorptionCurve = {{0, 0}, {1, 180}};
    const auto & solidMaterial = multiDomain.materials().createSolidMaterial(params);

    Gases::CGas gas;
    gas.addGasItem(0.1, Gases::GasDef::Air);
    gas.addGasItem(0.3, Gases::GasDef::Argon);
    gas.addGasItem(0.3, Gases::GasDef::Krypton);
    gas.addGasItem(0.3, Gases::GasDef::Xenon);

    auto & frameCavity = multiDomain.materials().createGas(
      "Frame Cavity 1", HygroThermFEM::CavityStandard::ISO15099, gas);

    // Elements that will contain frame cavity
    std::set<size_t> frameCavityElement{6, 7, 10};

    // Create elements grid
    size_t elementNumber{0u};
    for(auto ix = 1u; ix < gridX.size(); ++ix)
    {
        for(auto iy = 1u; iy < gridY.size(); ++iy)
        {
            {
                ++elementNumber;
                const auto node1 = ix * gridX.size() + iy - gridX.size();
                const auto node2 = ix * gridX.size() + iy - gridX.size() + 1u;
                const auto node3 = ix * gridX.size() + (iy + 1u);
                const auto node4 = ix * gridX.size() + iy;
                std::string materialName;
                if(frameCavityElement.find(elementNumber) != frameCavityElement.end())
                {
                    materialName = frameCavity.name();
                }
                else
                {
                    materialName = solidMaterial.name();
                }
                multiDomain.createElement({.node1 = node1, .node2 = node2, .node3 = node3, .node4 = node4, .material = materialName});
            }
        }
    }

    // Create Boundary Conditions
    const auto tAir{0.0};
    const auto hc{30.0};
    const auto humidity{0.6};

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff{tAir, hc, humidity};

    // Build boundary condition nodes on left edge
    std::vector<size_t> bcnodes;
    for(size_t i = 0u; i < gridY.size(); ++i)
    {
        bcnodes.push_back(i * gridX.size() + 1);
    }

    // Now build boundary condition on left edge of domain rectangle
    for(size_t i = 1u; i < bcnodes.size(); ++i)
    {
        multiDomain.createBC_FixedHc(bcnodes[i - 1u], bcnodes[i], bcCoeff);
    }

    // Now perform transient calculation in order to make frame cavity update over the simulation
    constexpr auto dTime = 36000;
    constexpr auto nSteps = 10;

    auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
    auto humidities = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);
    std::vector<std::vector<double>> temperatureSolution;
    std::vector<std::vector<double>> waterContentSolution;

    for(auto i = 0; i < nSteps; ++i)
    {
        auto aSolution = multiDomain.transient(temperatures, humidities, dTime);
        temperatureSolution.push_back(aSolution.temperature);
        waterContentSolution.push_back(aSolution.waterContent);
        temperatures = aSolution.temperature;
        humidities = aSolution.humidity;
    }

    // clang-format off
    std::vector<std::vector<double>> correctWaterContentSolution{{28.3368203561, 1.46351757768, 0.071560805145, 0.00559659926934, 0.00494633549068, 28.3381541817, 1.24268618438, 0.192539151833, 0.101234216742, 0.0129478776779, 28.3407365899, 1.03832086466, 0.275844991128, 0.15519774309, 0.0169595364437, 28.3381352886, 1.24043512741, 0.27529605839, 0.0204091929613, 0.00510007788118, 28.3368169862, 1.46323131201, 0.0700768467949, 0.00867657611898, 0.00118979890892},
 {71.6834672657, 2.85250619039, 0.122141014349, 0.00916649263333, 0.00789443809426, 72.3098585976, 2.47094411696, 0.338098334575, 0.167697843468, 0.0207263201747, 73.3861601109, 2.11422571075, 0.49042917165, 0.257990961016, 0.0273141946958, 71.9470733649, 2.46471159478, 0.476624965288, 0.0335127885211, 0.00820348455595, 71.3826526149, 2.85037900403, 0.120140353604, 0.014083088906, 0.00189187582562},
 {91.6179291363, 5.02829376948, 0.194960825952, 0.0141186698906, 0.0115414028985, 91.9778976097, 4.42503493732, 0.566083936471, 0.259984329646, 0.030265606194, 92.7269734783, 3.85035018619, 0.850567675496, 0.400390766233, 0.0398720841596, 91.7339352345, 4.39488222169, 0.78134382578, 0.0511066246347, 0.0119597789485, 91.425017454, 5.00336925476, 0.192020072254, 0.0211160995966, 0.00273771594093},
 {99.7151331166, 8.60575008085, 0.311932775474, 0.0219287363241, 0.0168958141005, 99.8925608817, 7.61177031147, 0.939462099806, 0.407308735867, 0.0442173853305, 100.348736209, 6.64209288053, 1.4470197734, 0.628155312421, 0.0581922357129, 99.7663406692, 7.55000865749, 1.27732614745, 0.0785671746123, 0.0174213130438, 99.6224109559, 8.55169243715, 0.307538033396, 0.0318066125877, 0.00394483393921},
 {102.825876087, 13.0461938264, 0.483213887892, 0.0337755913734, 0.0250289774617, 102.938836751, 11.4493850735, 1.46457628988, 0.631747612483, 0.0653744908606, 103.260376357, 9.87200074035, 2.26371739018, 0.975280662589, 0.0859621941926, 102.873090037, 11.3637404165, 1.98946189429, 0.120323786625, 0.0256810298724, 102.780328132, 12.9736532201, 0.47671657676, 0.0480271361364, 0.0057454817834},
 {104.032363221, 17.7164009279, 0.709513564102, 0.0503882134215, 0.037094851369, 104.138240966, 15.3370910409, 2.11139224604, 0.945352160331, 0.0967583956523, 104.424866599, 12.9943852038, 3.2219686726, 1.45969981704, 0.127159324831, 104.096797396, 15.2335069132, 2.8965135144, 0.179477435545, 0.0379276176675, 104.003679239, 17.6332503037, 0.699772165476, 0.0713987604026, 0.00840670736025},
 {104.542720947, 22.3057033214, 0.987362159648, 0.0720753428998, 0.0542227930144, 104.659829316, 19.0101219012, 2.84396264128, 1.35076865308, 0.141351239022, 104.955576838, 15.8105903887, 4.24217879653, 2.08435883518, 0.185703159681, 104.625749199, 18.8913942556, 3.96251682786, 0.257765099316, 0.0553476219411, 104.518346684, 22.2151286208, 0.972485175244, 0.103202586998, 0.0122164327259},
 {104.801156715, 26.7018057553, 1.31205951531, 0.0988591338567, 0.0774340289396, 104.931277869, 22.3943652612, 3.63277778957, 1.8440273632, 0.201880190977, 105.249840713, 18.2959869276, 5.26606674459, 2.84169577349, 0.265169492349, 104.897739064, 22.2628619335, 5.15210769199, 0.355861417956, 0.0790429296573, 104.776869775, 26.6049657249, 1.28932885026, 0.144391702846, 0.0174690011775},
 {104.967361212, 30.8737482131, 1.67885508732, 0.130606773121, 0.107597721812, 105.106860781, 25.4880519483, 4.45657860556, 2.4173666721, 0.280700220119, 105.448591954, 20.4816315076, 6.25728028186, 3.71817088898, 0.368634233034, 105.071580899, 25.3472083489, 6.43513893905, 0.473731915551, 0.109984313088, 104.942191569, 30.7713293814, 1.64477816306, 0.195643723975, 0.0244549995975},
 {105.098291633, 34.8198464047, 2.08333384847, 0.167123854199, 0.145416561248, 105.242446405, 28.3118330849, 5.30079807306, 3.06130817065, 0.379751480438, 105.602915115, 22.4091116314, 7.19589251014, 4.69760023118, 0.498616454363, 105.205008091, 28.1664889214, 7.78708164522, 0.610913278274, 0.148992936907, 105.072348757, 34.7125831774, 2.03376309149, 0.257416540471, 0.0334549524722}};
    // clang-format on

    // Inspecting only first ten nodes

    TestHelper::dumpGolden("correctWaterContentSolution", waterContentSolution);

    for(auto i = 0u; i < correctWaterContentSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctWaterContentSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctWaterContentSolution[i][j], waterContentSolution[i][j], 1e-6);
        }
    }
}
