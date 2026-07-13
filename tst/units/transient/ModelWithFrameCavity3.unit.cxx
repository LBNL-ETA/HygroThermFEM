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
 {71.6834655965, 2.85251703041, 0.122116515011, 0.00919773427034, 0.00788213911797, 72.3098359035, 2.47041295093, 0.339046684633, 0.167425321122, 0.0207011507991, 73.3861810977, 2.11325024007, 0.491489601089, 0.257758651151, 0.0272765374948, 71.9470375936, 2.4642112061, 0.477269669092, 0.0335035808402, 0.00818340007293, 71.3826601925, 2.85038130487, 0.12016383929, 0.0140800062229, 0.00189131923063},
 {91.617942653, 5.02831758267, 0.19491364128, 0.0141685461172, 0.0115245433082, 91.977881968, 4.42443009127, 0.567161083959, 0.259726536564, 0.0302286161537, 92.72695806, 3.84926043454, 0.85160730481, 0.400208751731, 0.0398202929663, 91.7339187267, 4.39428130816, 0.782136632807, 0.05110844089, 0.0119338859796, 91.425042713, 5.00337966879, 0.192041734272, 0.0211167146197, 0.00273730227833},
 {99.7151506324, 8.60579120707, 0.311865862281, 0.0219927814077, 0.0168766210277, 99.8925494604, 7.61117856172, 0.940537081591, 0.407101021388, 0.0441730478988, 100.348693878, 6.64103240408, 1.44791296077, 0.628045662168, 0.0581327530057, 99.7663349015, 7.54939153944, 1.27817856774, 0.0785808464906, 0.0173927077762, 99.6224389332, 8.55171913931, 0.307553542151, 0.0318116263359, 0.00394483491056},
 {102.825892548, 13.0462547338, 0.483129561111, 0.033852202294, 0.0250075394734, 102.938825748, 11.4487990859, 1.46568612736, 0.631563494981, 0.0653234669067, 103.260318798, 9.87093990011, 2.26454628786, 0.975207037328, 0.0858951855594, 102.873085779, 11.363104966, 1.99038885588, 0.120345204031, 0.025649197284, 102.780353029, 12.9737007167, 0.476726524673, 0.0480357830898, 0.00574599858467},
 {104.032378464, 17.7164827233, 0.709412657396, 0.0504785731762, 0.0370694610826, 104.138228519, 15.3364574764, 2.11265159864, 0.94514188602, 0.0966980679281, 104.42479837, 12.9932146485, 3.22289543243, 1.45960922382, 0.127080097449, 104.096791081, 15.2328021854, 2.89758599189, 0.179501336213, 0.0378894425966, 104.003701659, 17.6333197774, 0.699779653175, 0.0714097323796, 0.00840771693081},
 {104.542736642, 22.3058085014, 0.987242987849, 0.0721835651081, 0.0541900677894, 104.659814026, 19.0093762728, 2.84550473804, 1.35048195159, 0.141276073211, 104.955494737, 15.8091941681, 4.24336540725, 2.08420477826, 0.185602868079, 104.6257397, 18.8905567125, 3.96382196258, 0.257787891132, 0.0552977568508, 104.518369513, 22.2152214965, 0.972494526574, 0.1032143264, 0.0122178320872},
 {104.801174888, 26.7019392274, 1.31191765162, 0.098991826842, 0.0773892248354, 104.931258507, 22.3934568906, 3.6347110147, 1.8436284834, 0.201782521802, 105.249737895, 18.2942900964, 5.26762315441, 2.84145138601, 0.265036302056, 104.897726035, 22.2618422432, 5.15371785293, 0.355882773755, 0.0789743339775, 104.7768963, 26.6050851475, 1.28934474071, 0.144402486478, 0.0174706665162},
 {104.967383533, 30.8739168449, 1.67868392214, 0.130772289159, 0.107535145601, 105.106836463, 25.4869530927, 4.45896947797, 2.41683929087, 0.280571184566, 105.448462092, 20.4796112031, 6.25925187922, 3.71783101512, 0.368454611233, 105.071564288, 25.3459810236, 6.43709946767, 0.47375497243, 0.109888854829, 104.942224481, 30.7714799657, 1.64480528474, 0.195651648624, 0.0244568427336},
 {105.098319141, 34.8200580463, 2.08312543672, 0.167331524922, 0.145329971499, 105.24241668, 28.3105375486, 5.30367254982, 3.06065350859, 0.379581848558, 105.602754612, 22.4067871207, 7.19827144626, 4.69717706304, 0.498376285313, 105.204987979, 28.1650508125, 7.78941132126, 0.610944022619, 0.14886193239, 105.072389709, 34.7127701357, 2.03380612899, 0.257419487858, 0.0334569638688}};
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
