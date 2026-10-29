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
    std::vector<std::vector<double>> correctWaterContentSolution{{28.3368203561, 1.46351757768, 0.071560805145, 0.00559659926934, 0.00494633549068, 28.3381541817, 1.24268617869, 0.192539143662, 0.101234214678, 0.0129478776779, 28.3407365899, 1.0383209669, 0.27584475842, 0.155197731093, 0.0169595364437, 28.3381352886, 1.24043511301, 0.275296051158, 0.0204091929613, 0.00510007788118, 28.3368169862, 1.46323131201, 0.0700768467949, 0.00867657611898, 0.00118979890892},
 {71.7281056747, 2.85263872514, 0.122075087881, 0.00919807818308, 0.00787762098778, 72.315251641, 2.47013499061, 0.339258446282, 0.167355032248, 0.0206751691466, 73.268735555, 2.11168793237, 0.493069749375, 0.257514738435, 0.0272140796902, 71.9798891481, 2.46407269446, 0.477399799384, 0.0334802557864, 0.00816778466748, 71.4556129073, 2.85062370269, 0.120101378101, 0.0140715500469, 0.00188969172492},
 {91.648563086, 5.03106622684, 0.194878095344, 0.0141676225304, 0.0115179315208, 91.9885777898, 4.42415617882, 0.567169562598, 0.259659367427, 0.0301918172088, 92.650360369, 3.84284217232, 0.852781297374, 0.399951707756, 0.0397327696592, 91.7644593863, 4.39582654331, 0.782136335543, 0.0510794445305, 0.0119116518736, 91.4753004788, 5.00826604955, 0.19200461101, 0.0211037137278, 0.00273476276481},
 {99.7306681384, 8.61230369792, 0.311884373979, 0.0219894609829, 0.0168678202565, 99.9009162413, 7.61150765159, 0.940280487044, 0.4069771588, 0.0441270133871, 100.306161403, 6.62891121016, 1.44813688489, 0.627685228534, 0.0580259111659, 99.7854295511, 7.55379127617, 1.27805378413, 0.0785418210811, 0.0173652515808, 99.6480716096, 8.56316991388, 0.307630997266, 0.0317959484611, 0.00394153220165},
 {102.83399944, 13.0554607135, 0.483234277279, 0.0338468691374, 0.0249950779907, 102.944307368, 11.449843011, 1.46522736068, 0.631311261398, 0.0652640091732, 103.235393777, 9.85528666043, 2.26384152566, 0.974603557288, 0.0857630262814, 102.884531043, 11.3696988116, 1.99021532071, 0.12028863606, 0.0256149765813, 102.793926171, 12.9898034036, 0.476975739594, 0.048019186024, 0.00574186579364},
 {104.037951821, 17.727472841, 0.709621157299, 0.0504730186624, 0.0370509392462, 104.141830136, 15.3380011503, 2.11208205197, 0.944703064133, 0.0966170419279, 104.407029725, 12.9754938237, 3.22152488995, 1.45863075109, 0.126908605459, 104.104535002, 15.2407759971, 2.89744649376, 0.179420369317, 0.0378449049147, 104.01314399, 17.6525196123, 0.700232214272, 0.0713940048221, 0.00840256842009},
 {104.547921785, 22.3182359642, 0.987564679627, 0.0721805118224, 0.0541624814727, 104.662451901, 19.0112065317, 2.84489749831, 1.34981723798, 0.14116232456, 104.939755179, 15.7898634483, 4.24163456689, 2.08274558725, 0.185371549875, 104.632082667, 18.8994711378, 3.96378449523, 0.257677257086, 0.0552377443132, 104.527154425, 22.2369542386, 0.973169528529, 0.103201166572, 0.0122114032804},
 {104.806717026, 26.7157385602, 1.31235856012, 0.098994436862, 0.0773492770661, 104.933500369, 22.3954465098, 3.63410771652, 1.84271705538, 0.201622876977, 105.233886054, 18.2733764077, 5.26578345151, 2.83943897581, 0.264719913268, 104.903772335, 22.2715163557, 5.1538317443, 0.355739096356, 0.0788925373564, 104.78623024, 26.6292554755, 1.29025684428, 0.144393496836, 0.0174626514133},
 {104.97344539, 30.8890999756, 1.67924927259, 0.130783778838, 0.107479513634, 105.10897645, 25.4890451321, 4.45838414819, 2.41567539039, 0.280351398425, 105.431776055, 20.4570614414, 6.25747429529, 3.71522501549, 0.368024947907, 105.077765237, 25.3563467274, 6.43739853601, 0.4735766809, 0.109778282214, 104.95237389, 30.7981070968, 1.64596821958, 0.195648400996, 0.0244469310337},
 {105.104859623, 34.8366307215, 2.0838207717, 0.167354913533, 0.145255503362, 105.244584979, 28.3127156588, 5.30310084763, 3.05924207178, 0.379287380377, 105.585115856, 22.3825942548, 7.19665481363, 4.69396367477, 0.49780374446, 105.211482059, 28.1760748042, 7.78991868748, 0.610731040114, 0.148715316041, 105.083291647, 34.7418506201, 2.03523374969, 0.25742359651, 0.0334448729096}};
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
