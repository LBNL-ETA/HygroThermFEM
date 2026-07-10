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
    std::vector<std::vector<double>> correctWaterContentSolution{{28.3368204, 1.46351758, 0.0715608051, 0.00559659927, 0.00494633549, 28.3381542, 1.24268618, 0.19253915, 0.101234216, 0.0129478777, 28.3407366, 1.03832087, 0.275844968, 0.155197742, 0.0169595364, 28.3381353, 1.24043512, 0.275296057, 0.020409193, 0.00510007788, 28.336817, 1.46323131, 0.0700768468, 0.00867657612, 0.00118979891},
 {71.688495, 2.85253935, 0.122103084, 0.00919680145, 0.00788078768, 72.3103629, 2.47042147, 0.339048292, 0.167405415, 0.0206963782, 73.3764667, 2.11316808, 0.491594962, 0.25772033, 0.027268248, 71.9495136, 2.46422929, 0.477253645, 0.0334983611, 0.0081811972, 71.3891384, 2.8504106, 0.120149585, 0.0140778019, 0.00189096123},
 {91.63, 5.02927078, 0.194873081, 0.0141650001, 0.0115196673, 91.9793854, 4.42438767, 0.567122177, 0.259638952, 0.0302117328, 92.7027441, 3.84752564, 0.851837783, 0.400039776, 0.0397911151, 91.740384, 4.39464577, 0.782041212, 0.0510880136, 0.0119260827, 91.4411432, 5.0047629, 0.192004573, 0.0211088659, 0.00273605137},
 {99.7369883, 8.61132553, 0.311801549, 0.0219840461, 0.0168635314, 99.8958246, 7.61094235, 0.940371413, 0.406827074, 0.0441270995, 100.304548, 6.63194201, 1.44825637, 0.627500795, 0.0580521462, 99.7787295, 7.55152268, 1.27790893, 0.0785220181, 0.0173712875, 99.6518593, 8.55989855, 0.307523236, 0.0317919862, 0.00394162286},
 {102.862359, 13.0625132, 0.483074768, 0.033835183, 0.0249760734, 102.944668, 11.4482644, 1.46532842, 0.630866008, 0.0652096303, 103.185817, 9.84541368, 2.26528984, 0.973767643, 0.085690743, 102.894282, 11.3694746, 1.98989053, 0.120200886, 0.0255954461, 102.829607, 12.997877, 0.476796308, 0.0479932446, 0.00573857344},
 {104.086268, 17.7500282, 0.709439743, 0.0504509118, 0.0370020378, 104.1474, 15.3357592, 2.11206016, 0.943698271, 0.0964474323, 104.314688, 12.9423318, 3.2245519, 1.45654758, 0.126622309, 104.128679, 15.2461093, 2.8968857, 0.179202941, 0.037770084, 104.076086, 17.6831257, 0.700117234, 0.0713281768, 0.00839218596},
 {104.613543, 22.3620124, 0.987467325, 0.0721448197, 0.054061696, 104.672482, 19.0090197, 2.84461883, 1.3479645, 0.140789213, 104.811702, 15.7269455, 4.2463104, 2.07876876, 0.184705701, 104.668141, 18.9130605, 3.96302514, 0.257257656, 0.0550651268, 104.612406, 22.2981097, 0.973338352, 0.103074894, 0.0121885851},
 {104.886548, 26.7843168, 1.31249252, 0.0989439671, 0.0771689232, 104.947026, 22.3940652, 3.63342993, 1.83976375, 0.20093454, 105.078081, 18.1779927, 5.27187048, 2.83300737, 0.263468208, 104.949061, 22.2949478, 5.15299058, 0.355045082, 0.0785692497, 104.888609, 26.7254616, 1.29098635, 0.144186537, 0.0174205944},
 {105.064497, 30.9843104, 1.67979142, 0.130719676, 0.10718793, 105.124997, 25.4891323, 4.45717064, 2.41143101, 0.279218933, 105.255612, 20.3292913, 6.26446132, 3.70592095, 0.365953175, 105.12975, 25.3902858, 6.43665895, 0.472542829, 0.109244428, 105.067188, 30.9311126, 1.64756764, 0.195343307, 0.0243777028},
 {105.20458, 34.9590052, 2.08496637, 0.167280676, 0.14481885, 105.262213, 28.3147731, 5.30124103, 3.05358184, 0.377571237, 105.394721, 22.2242904, 7.20387081, 4.68151686, 0.494662412, 105.268028, 28.2205003, 7.78951869, 0.609301645, 0.147907259, 105.206741, 34.9114736, 2.03802361, 0.257007438, 0.0333397873}};
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
