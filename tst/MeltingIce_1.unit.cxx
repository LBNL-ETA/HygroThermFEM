#include <memory>
#include <gtest/gtest.h>

#include "MoisThermFEM2D.hxx"

using MoisThermFEM::NodePool;
using MoisThermFEM::MaterialPool;

class MeltingIce_1 : public testing::Test
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

TEST_F(MeltingIce_1, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Simple two elements example with moisture and heat transfer.");

    // Enter nodes. Arguments are: node number, x-coordinate, y-coordinate

    std::vector<double> gridXCoordinates{
      0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1};

    const double initialTemperature = -1.0;
    const double initialMoistureContent = 1.0;
    const double initialPressure = 101325;

    auto state =
      MoisThermFEM::State(initialTemperature, initialMoistureContent, initialPressure, 0);
    size_t nodeIndex = 0;
    for(auto val : gridXCoordinates)
    {
        ++nodeIndex;
        NodePool::Instance().createNode(nodeIndex, val, 0.00, state);
        ++nodeIndex;
        NodePool::Instance().createNode(nodeIndex, val, 0.05, state);
    }

    // Test material intentionally have lower density and specific heat so we can observe melting
    // ice
    auto & material = MaterialPool::Instance().createMaterial(
      "Test material",
      2,                          /// density
      0.22,                       /// porosity
      8,                          /// specific heat capacity (dry)
      15,                         /// diffusion resistance factor
      {{0.0, 1.8}, {180, 1.8}},   /// thermal conductivity as function of water content
      {{0, 0},                    /// liquid transportation coefficient
       {27, 1E-8},
       {45, 1.1E-8},
       {90, 2E-8},
       {126, 3.5E-8},
       {144, 5E-8},
       {162, 1E-7},
       {171, 2E-7},
       {180, 7E-7}},
      {{0, 0},   /// sorption curve
       {0.5, 5.3},
       {0.65, 8.4},
       {0.8, 12},
       {0.93, 17},
       {0.95, 25},
       {0.99, 63},
       {0.995, 83},
       {0.999, 120},
       {1, 180}});

    MoisThermFEM::MultiDomain domain;

    /// Create elements
    for(size_t i = 1; i <= (NodePool::Instance().maxIndex() - 2) / 2; ++i)
    {
        const auto node1 = 2u * i + 1u;
        const auto node2 = 2u * i + 2u;
        const auto node3 = 2u * i;
        const auto node4 = 2u * i - 1u;
        domain.createElement(node1, node2, node3, node4, material.name());
    }

    /// Create Boundary Conditions
    const auto hc = 5.0;
    const auto airTemperature = 20.0;
    const auto humidity = 0.0;

    domain.createConvectionBC(21, 22, hc, airTemperature, humidity);

    const auto dTime = 360;
    const auto nSteps = 10;

    auto temperatures = NodePool::Instance().properties(MoisThermFEM::Variable::temperature);
    auto humidities = NodePool::Instance().properties(MoisThermFEM::Variable::humidity);
    std::vector<std::vector<double>> temperatureSolution;
    std::vector<std::vector<double>> waterContentSolution;

    for(auto i = 0; i < nSteps; ++i)
    {
        auto aSolution = domain.transient(temperatures, humidities, dTime);
        temperatureSolution.push_back(aSolution.temperature);
        waterContentSolution.push_back(aSolution.waterContent);
        temperatures = aSolution.temperature;
        humidities = aSolution.humidity;
    }

    std::vector<std::vector<double>> correctWaterContentSolution{
      {180, 180, 180, 180, 180, 180, 180, 180, 180, 180, 178.9059595},
      {180, 180, 180, 180, 180, 180, 180, 180, 180, 180, 177.7804451},
      {180, 180, 180, 180, 180, 180, 180, 180, 180, 180, 176.6323806},
      {180, 180, 180, 180, 180, 180, 180, 180, 180, 180, 175.4659424},
      {180, 180, 180, 180, 180, 180, 180, 180, 180, 180, 174.2836191},
      {180, 180, 180, 180, 180, 180, 180, 180, 180, 180, 173.0870689},
      {180, 180, 180, 180, 180, 180, 180, 180, 180, 180, 171.8774618},
      {180, 180, 180, 180, 180, 180, 180, 180, 180, 180, 170.6556516},
      {180, 180, 180, 180, 180, 180, 180, 180, 180, 180, 169.4222743},
      {180, 180, 180, 180, 180, 180, 180, 180, 180, 180, 168.1778122}};

    for(size_t i = 0u; i < correctWaterContentSolution.size(); ++i)
    {
        for(size_t j = 0u; j < correctWaterContentSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctWaterContentSolution[i][j], waterContentSolution[i][2 * j], 1e-6);
        }
    }

    std::vector<std::vector<double>> correctTemperatureSolution{
        {-0.99058987, -0.98921071, -0.98466896, -0.97563333, -0.95945526, -0.93139259, -0.88321951, -0.80081537, -0.66002562, -0.41958155, -0.0090049373 },
        {-0.96121112, -0.95690533, -0.94313008, -0.91717883, -0.87409322, -0.80598604, -0.70111929, -0.54287478, -0.30902193, 0.027718199, 0.49556982 },
        {-0.90658904, -0.89858354, -0.87348259, -0.82796634, -0.75629986, -0.65010542, -0.49821878, -0.28685733, -0.00045124749, 0.37640401, 0.85546428 },
        {-0.82684747, -0.81516042, -0.77902014, -0.71519073, -0.61830418, -0.48096792, -0.29405356, -0.047293669, 0.26968786, 0.66585339, 1.14686 },
        {-0.72519324, -0.71029464, -0.66465746, -0.58549794, -0.46832248, -0.30718392, -0.095105258, 0.17528974, 0.51092907, 0.91728176, 1.3973303 },
        {-0.60586662, -0.58837793, -0.53515264, -0.4439665, -0.31129424, -0.13259337, 0.097284013, 0.38355513, 0.73087363, 1.142663, 1.6205128 },
        {-0.47298447, -0.45350906, -0.39450047, -0.29426349, -0.1501451, 0.041209814, 0.28351045, 0.58039835, 0.93498554, 1.3494026, 1.8244156 },
        {-0.33006703, -0.30912083, -0.24585111, -0.13900884, 0.013342154, 0.21361503, 0.46442386, 0.76826257, 1.1271686, 1.542408, 2.014217 },
        {-0.17993661, -0.15793326, -0.091613355, 0.019917197, 0.17803264, 0.38442265, 0.64088032, 0.94906145, 1.3102389, 1.7250785, 2.1934579 },
        {-0.024777176, -0.0020367721, 0.066400489, 0.18115525, 0.34317261, 0.5535963, 0.81360868, 1.1242517, 1.4862471, 1.8998345, 2.3646415 }};

    for(size_t i = 0u; i < correctTemperatureSolution.size(); ++i)
    {
        for(size_t j = 0u; j < correctTemperatureSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctTemperatureSolution[i][j], temperatureSolution[i][2 * j], 1e-6);
        }
    }
}
