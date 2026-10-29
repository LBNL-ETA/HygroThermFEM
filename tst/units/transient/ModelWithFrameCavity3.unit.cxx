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
    // Re-captured when the element started reading the cavity's conductivity live:
    // the moisture field feels the cavity's keff from then on.
    std::vector<std::vector<double>> correctWaterContentSolution{{28.3368203561, 1.46351757753, 0.0715608051422, 0.00559659926904, 0.00494633549029, 28.3381541817, 1.24268617858, 0.192539143646, 0.101234214669, 0.0129478776769, 28.3407365899, 1.03832096682, 0.275844758397, 0.15519773108, 0.0169595364425, 28.3381352887, 1.24043511289, 0.275296051135, 0.02040919296, 0.00510007788082, 28.3368169862, 1.46323131186, 0.0700768467922, 0.00867657611841, 0.00118979890885},
 {71.7281084573, 2.85263271124, 0.122096773146, 0.00917030963761, 0.00788774220412, 72.3152697553, 2.4705879077, 0.338465589558, 0.167576370154, 0.0206957482811, 73.2687135493, 2.11250814261, 0.492190086396, 0.257703948804, 0.0272444749478, 71.9799182015, 2.46449350107, 0.476868944491, 0.0334853002654, 0.00818379068629, 71.4556074737, 2.85062478824, 0.120079758977, 0.0140742422002, 0.00189010440729},
 {91.6485620652, 5.03102355637, 0.19491856566, 0.0141238171069, 0.011531667335, 91.9886004308, 4.42464610738, 0.566266592401, 0.259866007554, 0.0302219146691, 92.650379218, 3.84373306441, 0.851915742446, 0.400096806552, 0.0397744311511, 91.764482181, 4.39630780216, 0.781480997138, 0.0510742988782, 0.0119321572167, 91.4752894072, 5.00823733997, 0.191983482346, 0.0211032785842, 0.00273501701516},
 {99.7306611982, 8.61224023155, 0.31194089836, 0.0219335854293, 0.0168833488078, 99.9009323843, 7.61198361779, 0.939376220668, 0.407140671834, 0.0441629367487, 100.306198261, 6.62977805741, 1.44738940505, 0.627768561878, 0.0580736197141, 99.785441037, 7.55428326194, 1.27734601796, 0.0785262194538, 0.0173878122068, 99.6480565875, 8.56312205164, 0.307614123904, 0.0317918286714, 0.00394140629609},
 {102.833987016, 13.0553888899, 0.483304936532, 0.0337803063409, 0.0250124344732, 102.944317252, 11.4503252599, 1.46428851996, 0.631455262878, 0.0653053188049, 103.235437754, 9.85616786879, 2.26314135421, 0.974656831915, 0.0858168496437, 102.88453579, 11.3702174383, 1.9894406821, 0.120266324936, 0.0256401795267, 102.793907788, 12.9897465634, 0.476962676157, 0.0480121143281, 0.00574128328192},
 {104.037951238, 17.7273561895, 0.709704760997, 0.0503946573277, 0.0370717469922, 104.14185023, 15.3385066497, 2.11100713373, 0.944870229935, 0.0966661207952, 104.407089697, 12.9764612609, 3.22072978871, 1.45869897289, 0.126972769166, 104.104550539, 15.2413367238, 2.89654055811, 0.179395880858, 0.0378755753311, 104.013138513, 17.652417442, 0.700219673341, 0.0713852029307, 0.00840156011788},
 {104.547917935, 22.3180914888, 0.987662489415, 0.0720866369794, 0.0541898558529, 104.662472022, 19.0118078006, 2.84357169877, 1.35004996052, 0.1412240583, 104.939825014, 15.7910305532, 4.24060678488, 2.08286755814, 0.185453826985, 104.632098583, 18.9001457755, 3.96267146704, 0.257653509628, 0.0552786373745, 104.527145716, 22.2368247509, 0.973153320509, 0.1031921093, 0.0122100682987},
 {104.806707357, 26.7155695286, 1.31247434258, 0.0988791194374, 0.0773875592355, 104.933520858, 22.3961929273, 3.63243926743, 1.8430447994, 0.201703971901, 105.233970878, 18.274810508, 5.26443207425, 2.83963638469, 0.26483062038, 104.903788293, 22.272353011, 5.15245020388, 0.355716079849, 0.0789498755083, 104.786214789, 26.6291035976, 1.29023249639, 0.144385848554, 0.0174611061229},
 {104.973419865, 30.8889296268, 1.67938911405, 0.130639547286, 0.107533919534, 105.108990545, 25.4899838327, 4.45631945284, 2.41611104242, 0.280459585796, 105.431875078, 20.4587977857, 6.25576676107, 3.71550136263, 0.368175847912, 105.077773965, 25.3573900352, 6.43571280955, 0.473551389844, 0.109859241683, 104.952341053, 30.7979578301, 1.64593157438, 0.195644040547, 0.0244452612439},
 {105.104829332, 34.8364405554, 2.08399176932, 0.167173408783, 0.145331776976, 105.244603421, 28.3138395393, 5.30062026851, 3.05978357647, 0.379430708607, 105.585240197, 22.3846030104, 7.19460180243, 4.69430787422, 0.498007098866, 105.211493783, 28.1773148053, 7.78791521324, 0.610698021827, 0.148827564606, 105.083251936, 34.741688103, 2.03518042482, 0.257424633564, 0.0334430977662}};
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
