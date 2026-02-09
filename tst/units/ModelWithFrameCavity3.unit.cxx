#include <gtest/gtest.h>
#include <memory>
#include <set>

#include "HygroThermFEM2D.hxx"
#include "TestMaterials.hxx"

using HygroThermFEM::State;

TEST(TestModelWithFrameCavity3, TestSingleFrameCavity)
{
    SCOPED_TRACE("Begin Test: Model with single frame cavity.");

    // Create multi-domain before materials
    HygroThermFEM::MultiDomain multiDomain;

    std::vector<double> gridX{0, 0.01, 0.02, 0.03, 0.04};
    std::vector<double> gridY{0, 0.05, 0.1, 0.15, 0.2};

    const auto initialTemperature{20.0};
    const auto initialHumidity{0.0};
    const auto initialPressure{101325.0};

    const State state({
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
    std::vector<std::vector<double>> correctWaterContentSolution{
        {77.938809,1.414799,0.026126,0.001246,0.000456,78.751416,1.309171,0.101852,0.025597,0.001186,80.111687},
        {99.542013,2.967956,0.072601,0.003953,0.001683,99.847542,2.692281,0.264599,0.079626,0.004346,100.52708},
        {105.124529,4.528578,0.139676,0.008458,0.004191,105.200454,4.022253,0.474263,0.16836,0.01081,105.479106},
        {106.541084,6.057355,0.226191,0.014892,0.008485,106.549658,5.272354,0.717101,0.293288,0.021906,106.657446},
        {106.905733,7.545449,0.330737,0.023269,0.015024,106.900068,6.444155,0.983222,0.453124,0.038858,106.938998},
        {107.006832,8.992194,0.451963,0.033553,0.024205,107.000498,7.544902,1.266047,0.645222,0.062738,107.012181},
        {107.042011,10.398893,0.588637,0.045696,0.036358,107.037073,8.581978,1.561228,0.866398,0.094456,107.038297},
        {107.060612,11.767229,0.739635,0.059658,0.051749,107.056502,9.561875,1.865812,1.11338,0.134764,107.054104},
        {107.074882,13.098878,0.903931,0.075407,0.070586,107.070778,10.49018,2.177697,1.383028,0.184265,107.067587},
        {107.087824,14.395421,1.080583,0.092931,0.093025,107.083117,11.371734,2.495302,1.672447,0.24343,107.08032}};
    // clang-format on

    // Inspecting only first ten nodes

    for(auto i = 0u; i < correctWaterContentSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctWaterContentSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctWaterContentSolution[i][j], waterContentSolution[i][j], 1e-6);
        }
    }
}
