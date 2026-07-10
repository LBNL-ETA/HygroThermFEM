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
    std::vector<std::vector<double>> correctWaterContentSolution{{28.3368203561, 1.46351757768, 0.071560805145, 0.00559659926934, 0.00494633549068, 28.3381541817, 1.24268618391, 0.192539151524, 0.10123421668, 0.0129478776779, 28.3407365899, 1.03832086425, 0.275844988484, 0.155197742929, 0.0169595364437, 28.3381352886, 1.24043512685, 0.275296058213, 0.0204091929613, 0.00510007788118, 28.3368169862, 1.46323131201, 0.0700768467949, 0.00867657611898, 0.00118979890892},
 {71.6843576277, 2.85252160021, 0.122114720679, 0.00919759578485, 0.00788195488367, 72.3103835709, 2.47041838223, 0.339045812436, 0.167422328249, 0.0207005425731, 73.3860428223, 2.11324961176, 0.491498371842, 0.257753203551, 0.0272755270405, 71.9478033245, 2.46421770434, 0.477266092177, 0.0335028501098, 0.00818312358715, 71.3837371957, 2.8503865814, 0.12016193308, 0.0140796992735, 0.00189127063106},
 {91.6189706413, 5.02841882682, 0.194911185656, 0.0141682874652, 0.0115241588917, 91.9785096052, 4.4244951387, 0.567161216046, 0.259721009078, 0.0302273616146, 92.7266779585, 3.84926226259, 0.851625836808, 0.400198537302, 0.0398182195409, 91.7348209195, 4.39436859874, 0.782132921236, 0.0511070378496, 0.0119333135873, 91.4263146034, 5.00350543307, 0.192039422547, 0.021116104845, 0.00273719917101},
 {99.7159864551, 8.60609945147, 0.311865409582, 0.0219924985227, 0.0168760449643, 99.8930624799, 7.61135958911, 0.940543169566, 0.407094248303, 0.0441711476687, 100.348388648, 6.64102278892, 1.44794526916, 0.628032357166, 0.0581295822135, 99.7670856821, 7.54964406894, 1.27818331522, 0.0785790051978, 0.0173918318092, 99.6234847672, 8.55211033773, 0.307554291409, 0.0318108218237, 0.00394467743158},
 {102.826544026, 13.0467799887, 0.483133710251, 0.0338520292124, 0.0250067824486, 102.939238538, 11.4490978507, 1.46570250269, 0.631556880407, 0.0653208955835, 103.260085026, 9.87091656492, 2.26459586607, 0.975192202155, 0.0858908001387, 102.873695264, 11.3635286834, 1.99041024273, 0.12034320323, 0.0256479924258, 102.78117578, 12.9743748766, 0.476733771669, 0.0480349429992, 0.0057457877958},
 {104.032927089, 17.717199541, 0.709423509167, 0.0504786619284, 0.037068534939, 104.138609571, 15.3368697709, 2.11268174366, 0.945136927233, 0.096694776011, 104.424733991, 12.9932142639, 3.22296554228, 1.45959453181, 0.127074312294, 104.097342682, 15.2333846485, 2.89763093793, 0.17949948614, 0.0378878695157, 104.004395752, 17.6342436731, 0.699796191364, 0.0714090473876, 0.00840745538164},
 {104.543252257, 22.3067002591, 0.987262279418, 0.0721840791066, 0.0541889849403, 104.660219506, 19.0099113986, 2.84555240731, 1.35048037798, 0.141271983465, 104.955665006, 15.8092754099, 4.24346146172, 2.0841923754, 0.185595428851, 104.626302233, 18.8912981329, 3.96389733372, 0.257786540783, 0.055295762328, 104.519014842, 22.2163687295, 0.972522627179, 0.103214008097, 0.0122175239007},
 {104.801697691, 26.7029999541, 1.31194695821, 0.0989929419075, 0.0773880067215, 104.931718108, 22.3941304816, 3.63478047575, 1.84363262701, 0.201777560754, 105.250170249, 18.2945116699, 5.26775211104, 2.84144467999, 0.265026952114, 104.898334273, 22.2627500792, 5.15383107308, 0.355882368104, 0.0789718684505, 104.777535738, 26.6064401524, 1.28938637742, 0.144402774217, 0.0174703190547},
 {104.967928826, 30.8751435716, 1.67872476739, 0.130774192739, 0.107533838747, 105.107358568, 25.4877783793, 4.45906510859, 2.41685232943, 0.280565343502, 105.449156202, 20.4800202847, 6.25942022144, 3.71783519474, 0.368443212018, 105.072225505, 25.3470603635, 6.43725830373, 0.473756106412, 0.109885903837, 104.94287063, 30.7730295009, 1.64486224938, 0.195652824144, 0.0244564699692},
 {105.098887645, 34.8214458088, 2.08317931611, 0.167334410927, 0.145328664019, 105.242997516, 28.3115213175, 5.30379840363, 3.06067941552, 0.37957524592, 105.603693959, 22.407416082, 7.19848377336, 4.69719899984, 0.498362933777, 105.20569444, 28.1662999487, 7.78962334398, 0.610947455037, 0.148858547595, 105.07303873, 34.7144980133, 2.03388004249, 0.257421886102, 0.0334565902173}};
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
