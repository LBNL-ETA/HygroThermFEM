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
    std::vector<std::vector<double>> correctWaterContentSolution{{28.3368203561, 1.46351757753, 0.0715608051422, 0.00559659926904, 0.00494633549029, 28.3381541817, 1.24268618427, 0.192539151817, 0.101234216734, 0.0129478776769, 28.3407365899, 1.03832086458, 0.275844991105, 0.155197743077, 0.0169595364425, 28.3381352887, 1.24043512729, 0.275296058367, 0.02040919296, 0.00510007788082, 28.3368169862, 1.46323131186, 0.0700768467922, 0.00867657611841, 0.00118979890885},
 {71.6834672682, 2.85250617405, 0.122141013695, 0.00916649258551, 0.00789443806636, 72.3098586092, 2.4709440757, 0.33809833008, 0.167697841998, 0.0207263200989, 73.3861601455, 2.1142256268, 0.490429160797, 0.257990958499, 0.0273141945896, 71.947073375, 2.46471155689, 0.476624959761, 0.0335127882821, 0.00820348452565, 71.3826526174, 2.85037898873, 0.120140352996, 0.0140830888486, 0.00189187582108},
 {91.6179377143, 5.02826727543, 0.194960384244, 0.0141186487231, 0.0115413959471, 91.9779075008, 4.42500724136, 0.566081912897, 0.259983876801, 0.0302655882789, 92.7269856445, 3.85032106824, 0.850563862749, 0.400390047422, 0.0398720603899, 91.7339440602, 4.39485706234, 0.781341538209, 0.0511065538815, 0.0119597721301, 91.4250250958, 5.00334560529, 0.192019670843, 0.0211160791665, 0.00273771479548},
 {99.7151392565, 8.60571671099, 0.311931854886, 0.0219286821475, 0.016895790891, 99.8925671925, 7.61173791529, 0.939458417379, 0.407307584686, 0.0442173253759, 100.348743129, 6.64206153101, 1.44701323155, 0.628153490741, 0.0581921562608, 99.766346449, 7.54997905325, 1.27732173178, 0.0785669889359, 0.0174212901818, 99.6224166519, 8.55166231321, 0.30753719594, 0.031806555511, 0.00394482998192},
 {102.825876075, 13.0461681311, 0.48321264437, 0.0337755037685, 0.0250289292204, 102.938836803, 11.4493610439, 1.46457194093, 0.631745752167, 0.0653743658037, 103.260376541, 9.87197831539, 2.26371017248, 0.975277723618, 0.085962028465, 102.873090048, 11.3637187909, 1.98945632603, 0.120323475573, 0.0256809820395, 102.780328122, 12.9736304047, 0.47671545724, 0.0480270346686, 0.00574547327206},
 {104.032374127, 17.7163475244, 0.709511637062, 0.0503880789835, 0.0370947665246, 104.138249815, 15.337045559, 2.1113860384, 0.94534935969, 0.0967581751862, 104.424873004, 12.9943477886, 3.22195890827, 1.45969541393, 0.127159032886, 104.096806295, 15.2334636677, 2.89650538983, 0.179476953078, 0.0379275331274, 104.003690219, 17.6331996181, 0.6997704017, 0.0713985945019, 0.00840669188222},
 {104.542728978, 22.3056413541, 0.987359468735, 0.072075149885, 0.0542226573275, 104.659835655, 19.0100715826, 2.84395469903, 1.35076473281, 0.14135088576, 104.955580882, 15.8105515203, 4.24216704253, 2.08435270797, 0.18570269232, 104.625755626, 18.8913460369, 3.96250612233, 0.25776440331, 0.0553474861032, 104.518354798, 22.2150692623, 0.972482689925, 0.103202335191, 0.0122164070907},
 {104.801160874, 26.7017426678, 1.31205608853, 0.0988588760852, 0.0774338282631, 104.931280973, 22.3943161337, 3.6327685418, 1.84402226628, 0.201879667578, 105.249842139, 18.2959508954, 5.2660539656, 2.84168785608, 0.265168800634, 104.897742271, 22.2628147946, 5.152094802, 0.355860483608, 0.0790427278243, 104.776874013, 26.6049051748, 1.28932567796, 0.144391347812, 0.0174689618997},
 {104.967352768, 30.8737129795, 1.67885134488, 0.130606462357, 0.107597446664, 105.106853083, 25.4880280059, 4.45656981562, 2.41736068671, 0.280699500874, 105.448584987, 20.4816163746, 6.25726954809, 3.71816164793, 0.368633283423, 105.071573449, 25.3471859166, 6.43512589343, 0.473730772712, 0.109984034833, 104.942183373, 30.7712961512, 1.64477472937, 0.195643265924, 0.0244549436888},
 {105.09828238, 34.8198272251, 2.08333002466, 0.167123506389, 0.145416208449, 105.242438216, 28.3118227987, 5.30079031614, 3.06130165922, 0.379750555737, 105.602908289, 22.4091066584, 7.1958844785, 4.69759024258, 0.498615234756, 105.205000139, 28.1664797475, 7.78706934005, 0.610911974638, 0.148992577926, 105.072339732, 34.7125656073, 2.03375963456, 0.257415988801, 0.033454877919}};
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
