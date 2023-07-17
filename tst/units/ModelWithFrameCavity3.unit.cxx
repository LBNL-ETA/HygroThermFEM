#include <gtest/gtest.h>
#include <memory>
#include <set>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::NodePool;
using HygroThermFEM::MaterialPool;
using HygroThermFEM::State;

class TestModelWithFrameCavity3 : public testing::Test
{
protected:
    void SetUp() override
    {}

    void TearDown() override
    {
        NodePool::Instance().clear();
        MaterialPool::Instance().clear();
    }
};

TEST_F(TestModelWithFrameCavity3, TestSingleFrameCavity)
{
    SCOPED_TRACE("Begin Test: Model with single frame cavity.");

    std::vector<double> gridX{0, 0.01, 0.02, 0.03, 0.04};
    std::vector<double> gridY{0, 0.05, 0.1, 0.15, 0.2};

    const auto initialTemperature{20.0};
    const auto initialHumidity{0.0};
    const auto initialPressure{101325.0};

    const State state(initialTemperature, initialHumidity, initialPressure);
    size_t nodeIndex = 0;

    // Crating grid nodes
    for(auto y : gridY)
    {
        for(auto x : gridX)
        {
            ++nodeIndex;
            NodePool::Instance().createNode(nodeIndex, x, y, state);
        }
    }

    // Material Properties
    const double thermalConductivityDry{1.8};
    const double density{2050.0};
    const double porosity{0.22};
    const double specificHeatCapacityDry{850.0};
    const double diffusionResistanceFactor{15.0};
    const std::vector<FenestrationCommon::point> thermalConductivityMoistureDependent = {
      {0.0, 1.8}, {180, 1.8}};
    const double thermalConductivityMeasuredAtTemperature{0};
    const std::vector<FenestrationCommon::point> thermalConductivityTemperatureDependent = {
      {0.0, 1.8}, {1, 1.8}};
    const double thermalConductivityMeasuredAtHumidity{0};
    const std::vector<FenestrationCommon::point> liquidTransportationCurve = {{0, 0}, {180, 2e-6}};

    const std::vector<FenestrationCommon::point> moistureStorageFunction = {{0, 0}, {1, 180}};

    auto & solidMaterial =
      MaterialPool::Instance().createSolidMaterial("Material 1",
                                                   thermalConductivityDry,
                                                   density,
                                                   porosity,
                                                   specificHeatCapacityDry,
                                                   diffusionResistanceFactor,
                                                   thermalConductivityMoistureDependent,
                                                   thermalConductivityMeasuredAtTemperature,
                                                   thermalConductivityTemperatureDependent,
                                                   thermalConductivityMeasuredAtHumidity,
                                                   liquidTransportationCurve,
                                                   moistureStorageFunction);

    Gases::CGas gas;
    gas.addGasItem(0.1, Gases::GasDef::Air);
    gas.addGasItem(0.3, Gases::GasDef::Argon);
    gas.addGasItem(0.3, Gases::GasDef::Krypton);
    gas.addGasItem(0.3, Gases::GasDef::Xenon);

    auto & frameCavity = MaterialPool::Instance().createGas(
      "Frame Cavity 1", HygroThermFEM::CavityStandard::ISO15099, gas);

    // Elements that will contain frame cavity
    std::set<size_t> frameCavityElement{6, 7, 10};

    // Create elements grid
    HygroThermFEM::MultiDomain domain;
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
                createElement(domain, node1, node2, node3, node4, materialName);
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
        createBC_FixedHc(domain, bcnodes[i - 1u], bcnodes[i], bcCoeff);
    }

    // Now perform transient calculation in order to make frame cavity update over the simulation
    const auto dTime = 36000;
    const auto nSteps = 10;

    auto temperatures = properties(HygroThermFEM::Variable::temperature);
    auto humidities = properties(HygroThermFEM::Variable::humidity);
    std::vector<std::vector<double>> temperatureSolution;
    std::vector<std::vector<double>> waterContentSolution;

    HygroThermFEM::TransientSubstitutionSolver solver;
    for(auto i = 0; i < nSteps; ++i)
    {
        auto aSolution = solver.transient(domain, temperatures, humidities, dTime, i);
        temperatureSolution.push_back(aSolution.temperature);
        waterContentSolution.push_back(aSolution.waterContent);
        temperatures = aSolution.temperature;
        humidities = aSolution.humidity;
    }

    // clang-format off
    std::vector<std::vector<double>> correctWaterContentSolution{
        {77.938728,1.414799,0.026126,0.001247,0.000456,78.751333,1.309171,0.101853,0.025598,0.001186,80.111600,1.184037,0.183176,0.040911,0.001607,78.278705,1.311282,0.126088,0.004245,0.000463,77.537547,1.416376,0.026539,0.001368,0.000083},
        {99.541842,2.967957,0.072601,0.003953,0.001683,99.847367,2.692281,0.264600,0.079626,0.004346,100.526897,2.385262,0.455709,0.125584,0.005788,99.623506,2.690245,0.336663,0.013459,0.001680,99.374416,2.968308,0.073342,0.004528,0.000313},
        {105.124331,4.528578,0.139677,0.008459,0.004191,105.200253,4.022253,0.474264,0.168361,0.010811,105.478898,3.484321,0.784083,0.263826,0.014297,105.123739,4.014954,0.622895,0.029146,0.004172,105.080275,4.527927,0.140452,0.010214,0.000799},
        {106.540858,6.057356,0.226192,0.014893,0.008485,106.549430,5.272353,0.717102,0.293289,0.021907,106.657215,4.472988,1.137971,0.457662,0.028883,106.527210,5.259615,0.971409,0.052121,0.008461,106.533814,6.055987,0.226379,0.018995,0.001664},
        {106.905474,7.545448,0.330738,0.023269,0.015025,106.899809,6.444153,0.983223,0.453125,0.038859,106.938742,5.366392,1.497473,0.704616,0.051150,106.895700,6.426474,1.370156,0.082753,0.015041,106.906953,7.543533,0.329350,0.031303,0.003035},
        {107.006506,8.992193,0.451964,0.033553,0.024206,107.000175,7.544899,1.266049,0.645224,0.062739,107.011877,6.180646,1.851115,1.000029,0.082505,107.001961,7.523504,1.809327,0.121148,0.024351,107.009077,8.989901,0.447682,0.047466,0.005040},
        {107.041841,10.398891,0.588638,0.045697,0.036359,107.036904,8.581974,1.561230,0.866400,0.094457,107.038141,6.928804,2.192894,1.338395,0.124134,107.040572,8.558740,2.281059,0.167262,0.036775,107.043895,10.396430,0.579857,0.067732,0.007804},
        {107.060358,11.767227,0.739636,0.059658,0.051750,107.056249,9.561871,1.865814,1.113382,0.134765,107.053862,7.620979,2.520082,1.714083,0.177012,107.060470,9.539191,2.778995,0.220969,0.052639,107.061666,11.764830,0.724531,0.092288,0.011449},
        {107.074665,13.098874,0.903933,0.075408,0.070587,107.070559,10.490175,2.177699,1.383030,0.184266,107.067363,8.265083,2.831826,2.121712,0.241917,107.074883,10.470815,3.297917,0.282101,0.072215,107.075247,13.096787,0.880514,0.121267,0.016093},
        {107.087636,14.395417,1.080585,0.092931,0.093026,107.082927,11.371727,2.495304,1.672450,0.243431,107.080111,8.867444,3.128300,2.556326,0.319446,107.087205,11.358693,3.833486,0.350467,0.095727,107.087542,14.393876,1.046751,0.154764,0.021847},
    };
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
