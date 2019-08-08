#include <memory>
#include <vector>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::NodePool;
using HygroThermFEM::MaterialPool;

class MultiDomain_2D_THERMM200mmSlab : public testing::Test
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

TEST_F(MultiDomain_2D_THERMM200mmSlab, TestExample_1)
{
    // Enter nodes. Arguments are: node number, x-coordinate, y-coordinate

    std::vector<double> gridXCoordinates{0, 0.05, 0.1};

    // const double initialTemperature = 0.0;
    // const double initialMoistureContent = 0.0;
    // const double initialPressure = 101325;
    // const auto liquidPercent = 1.0;

    // auto state = HygroThermFEM::State(
    //  initialTemperature, initialMoistureContent, initialPressure, liquidPercent);

    NodePool::Instance().createNode(1, 0.1, -0.049);
    NodePool::Instance().createNode(2, 0.1, -0.009);
    NodePool::Instance().createNode(3, 0.06, -0.049);
    NodePool::Instance().createNode(4, 0.1, 0.049);
    NodePool::Instance().createNode(5, 0.06, -0.009);
    NodePool::Instance().createNode(6, 0, -0.049);
    NodePool::Instance().createNode(7, 0.06, 0.049);
    NodePool::Instance().createNode(8, 0, -0.009);
    NodePool::Instance().createNode(9, -0.06, -0.049);
    NodePool::Instance().createNode(10, 0, 0.049);
    NodePool::Instance().createNode(11, -0.06, -0.009);
    NodePool::Instance().createNode(12, -0.1, -0.049);
    NodePool::Instance().createNode(13, -0.06, 0.049);
    NodePool::Instance().createNode(14, -0.1, -0.009);
    NodePool::Instance().createNode(15, -0.1, 0.049);

    auto & material = MaterialPool::Instance().createSolidMaterial(
      "Aerated Concrete (density: 400 kg/m)",
      400,                        /// density
      0.81,                       /// porosity
      850,                        /// specific heat capacity (dry)
      7.9,                        /// diffusion resistance factor
      {{0.0, 0.1}, {380, 0.1}},   /// thermal conductivity as function of water content
      {{0, 0},                    /// liquid transportation coefficient
       {3.1, 2e-10},
       {38, 5.4e-09},
       {265, 1.1e-08},
       {342, 2e-08},
       {380, 1e-07}},
      {{0, 0},
       {0.5, 3.1},
       {0.8, 8.4},
       {0.93, 18},
       {0.964, 116},
       {0.99, 266},
       {0.995, 269},
       {0.999, 277},
       {1, 380}});

    HygroThermFEM::MultiDomain domain;

    /// Create elements
    domain.createElement(3, 1, 2, 5, material.name());
    domain.createElement(2, 4, 7, 5, material.name());
    domain.createElement(6, 3, 5, 8, material.name());
    domain.createElement(5, 7, 10, 8, material.name());
    domain.createElement(9, 6, 8, 11, material.name());
    domain.createElement(11, 14, 12, 9, material.name());
    domain.createElement(11, 8, 10, 13, material.name());
    domain.createElement(11, 13, 15, 14, material.name());

    /// Create Boundary Conditions
    const auto hc = 2.0;
    const auto airTemperature = 20.0;
    const auto humidity = 0.4;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff{airTemperature, hc, humidity};
    const std::vector<HygroThermFEM::FixedBCHCCoefficients> bcCoeffs{
      bcCoeff, bcCoeff, bcCoeff, bcCoeff, bcCoeff, bcCoeff, bcCoeff, bcCoeff, bcCoeff, bcCoeff,
      bcCoeff, bcCoeff, bcCoeff, bcCoeff, bcCoeff, bcCoeff, bcCoeff, bcCoeff, bcCoeff, bcCoeff};

    domain.createMoistureBCFixedHc(15, 14, bcCoeffs);
    domain.createMoistureBCFixedHc(14, 12, bcCoeffs);
    domain.createMoistureBCFixedHc(1, 2, bcCoeffs);
    domain.createMoistureBCFixedHc(2, 4, bcCoeffs);

    const auto dTime = 3600;
    const auto nSteps = 20;

    auto temperatures = NodePool::Instance().properties(HygroThermFEM::Variable::temperature);
    auto humidities = NodePool::Instance().properties(HygroThermFEM::Variable::humidity);
    std::vector<std::vector<double>> temperatureSolution;
    std::vector<std::vector<double>> humiditySolution;
    size_t timestepIndex{0u};

    for(auto i = 0; i < nSteps; ++i)
    {
        auto aSolution = domain.transient(temperatures, humidities, dTime, timestepIndex);
        temperatureSolution.push_back(aSolution.temperature);
        humiditySolution.push_back(aSolution.humidity);
        temperatures = aSolution.temperature;
        humidities = aSolution.humidity;
        ++timestepIndex;
    }

    // clang-format off
    std::vector<std::vector<double>> correctHumiditySolution{
    {0.212371, 0.212371, 0.001157, 0.212371, 0.001157, 6e-06, 0.001157, 6e-06, 0.001157, 6e-06, 0.001157, 0.212371, 0.001157, 0.212371, 0.212371},
    {0.317107, 0.317107, 0.003191, 0.317107, 0.003191, 2.5e-05, 0.003191, 2.5e-05, 0.003191, 2.5e-05, 0.003191, 0.317107, 0.003191, 0.317107, 0.317107},
    {0.373775, 0.373775, 0.00589, 0.373775, 0.00589, 6.5e-05, 0.00589, 6.5e-05, 0.00589, 6.5e-05, 0.00589, 0.373775, 0.00589, 0.373775, 0.373775},
    {0.40718, 0.40718, 0.009101, 0.40718, 0.009101, 0.000134, 0.009101, 0.000134, 0.009101, 0.000134, 0.009101, 0.40718, 0.009101, 0.40718, 0.40718},
    {0.427134, 0.427134, 0.012703, 0.427134, 0.012703, 0.000242, 0.012703, 0.000242, 0.012703, 0.000242, 0.012703, 0.427134, 0.012703, 0.427134, 0.427134},
    {0.438536, 0.438536, 0.016595, 0.438536, 0.016595, 0.000396, 0.016595, 0.000396, 0.016595, 0.000396, 0.016595, 0.438536, 0.016595, 0.438536, 0.438536},
    {0.44423, 0.44423, 0.020698, 0.44423, 0.020698, 0.000603, 0.020698, 0.000603, 0.020698, 0.000603, 0.020698, 0.44423, 0.020698, 0.44423, 0.44423},
    {0.446048, 0.446048, 0.024947, 0.446048, 0.024947, 0.000868, 0.024947, 0.000868, 0.024947, 0.000868, 0.024947, 0.446048, 0.024947, 0.446048, 0.446048},
    {0.445239, 0.445239, 0.02929, 0.445239, 0.02929, 0.001195, 0.02929, 0.001195, 0.02929, 0.001195, 0.02929, 0.445239, 0.02929, 0.445239, 0.445239},
    {0.442689, 0.442689, 0.033685, 0.442689, 0.033685, 0.001588, 0.033685, 0.001588, 0.033685, 0.001588, 0.033685, 0.442689, 0.033685, 0.442689, 0.442689},
    {0.439035, 0.439035, 0.038098, 0.439035, 0.038098, 0.002049, 0.038098, 0.002049, 0.038098, 0.002049, 0.038098, 0.439035, 0.038098, 0.439035, 0.439035},
    {0.434744, 0.434744, 0.042504, 0.434744, 0.042504, 0.002577, 0.042504, 0.002577, 0.042504, 0.002577, 0.042504, 0.434744, 0.042504, 0.434744, 0.434744},
    {0.430149, 0.430149, 0.046882, 0.430149, 0.046882, 0.003175, 0.046882, 0.003175, 0.046882, 0.003175, 0.046882, 0.430149, 0.046882, 0.430149, 0.430149},
    {0.425488, 0.425488, 0.051217, 0.425488, 0.051217, 0.003841, 0.051217, 0.003841, 0.051217, 0.003841, 0.051217, 0.425488, 0.051217, 0.425488, 0.425488},
    {0.420927, 0.420927, 0.055497, 0.420927, 0.055497, 0.004574, 0.055497, 0.004574, 0.055497, 0.004574, 0.055497, 0.420927, 0.055497, 0.420927, 0.420927},
    {0.416572, 0.416572, 0.059714, 0.416572, 0.059714, 0.005374, 0.059714, 0.005374, 0.059714, 0.005374, 0.059714, 0.416572, 0.059714, 0.416572, 0.416572},
    {0.412491, 0.412491, 0.063863, 0.412491, 0.063863, 0.006239, 0.063863, 0.006239, 0.063863, 0.006239, 0.063863, 0.412491, 0.063863, 0.412491, 0.412491},
    {0.40872, 0.40872, 0.067939, 0.40872, 0.067939, 0.007166, 0.067939, 0.007166, 0.067939, 0.007166, 0.067939, 0.40872, 0.067939, 0.40872, 0.40872},
    {0.405274, 0.405274, 0.07194, 0.405274, 0.07194, 0.008155, 0.07194, 0.008155, 0.07194, 0.008155, 0.07194, 0.405274, 0.07194, 0.405274, 0.405274},
    {0.402152, 0.402152, 0.075864, 0.402152, 0.075864, 0.009202, 0.075864, 0.009202, 0.075864, 0.009202, 0.075864, 0.402152, 0.075864, 0.402152, 0.402152}};
    // clang-format on

    EXPECT_EQ(humiditySolution.size(), correctHumiditySolution.size());
    
    for(auto i = 0u; i < humiditySolution.size(); ++i)
    {
        for(auto j = 0u; j < humiditySolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctHumiditySolution[i][j], humiditySolution[i][j], 1e-6);
        }
    }

    // clang-format off
    std::vector<std::vector<double>> correctTemperatureSolution{
    {10.422203, 10.422203, 3.110194, 10.422203, 3.110194, 1.151756, 3.110194, 1.151756, 3.110194, 1.151756, 3.110194, 10.422203, 3.110194, 10.422203, 10.422203},
    {13.656361, 13.656361, 5.977402, 13.656361, 5.977402, 2.938471, 5.977402, 2.938471, 5.977402, 2.938471, 5.977402, 13.656361, 5.977402, 13.656361, 13.656361},
    {14.900766, 14.900766, 8.198621, 14.900766, 8.198621, 4.885485, 8.198621, 4.885485, 8.198621, 4.885485, 8.198621, 14.900766, 8.198621, 14.900766, 14.900766},
    {15.569699, 15.569699, 9.904937, 15.569699, 9.904937, 6.742509, 9.904937, 6.742509, 9.904937, 6.742509, 9.904937, 15.569699, 9.904937, 15.569699, 15.569699},
    {16.047299, 16.047299, 11.25223, 16.047299, 11.25223, 8.409717, 11.25223, 8.409717, 11.25223, 8.409717, 11.25223, 16.047299, 11.25223, 16.047299, 16.047299},
    {16.444314, 16.444314, 12.348742, 16.444314, 12.348742, 9.864378, 12.348742, 9.864378, 12.348742, 9.864378, 12.348742, 16.444314, 12.348742, 16.444314, 16.444314},
    {16.796557, 16.796557, 13.263386, 16.796557, 13.263386, 11.117729, 13.263386, 11.117729, 13.263386, 11.117729, 13.263386, 16.796557, 13.263386, 16.796557, 16.796557},
    {17.117624, 17.117624, 14.040542, 17.117624, 14.040542, 12.193311, 14.040542, 12.193311, 14.040542, 12.193311, 14.040542, 17.117624, 14.040542, 17.117624, 17.117624},
    {17.413448, 17.413448, 14.709981, 17.413448, 14.709981, 13.116998, 14.709981, 13.116998, 14.709981, 13.116998, 14.709981, 17.413448, 14.709981, 17.413448, 17.413448},
    {17.68708, 17.68708, 15.292654, 17.68708, 15.292654, 13.912861, 15.292654, 13.912861, 15.292654, 13.912861, 15.292654, 17.68708, 15.292654, 17.68708, 17.68708},
    {17.940307, 17.940307, 15.803926, 17.940307, 15.803926, 14.601767, 15.803926, 14.601767, 15.803926, 14.601767, 15.803926, 17.940307, 15.803926, 17.940307, 17.940307},
    {18.174341, 18.174341, 16.255438, 18.174341, 16.255438, 15.201184, 16.255438, 15.201184, 16.255438, 15.201184, 16.255438, 18.174341, 16.255438, 18.174341, 18.174341},
    {18.390137, 18.390137, 16.656228, 18.390137, 16.656228, 15.72547, 16.656228, 15.72547, 16.656228, 15.72547, 16.656228, 18.390137, 16.656228, 18.390137, 18.390137},
    {18.588554, 18.588554, 17.013448, 18.588554, 17.013448, 16.18633, 17.013448, 16.18633, 17.013448, 16.18633, 17.013448, 18.588554, 17.013448, 18.588554, 18.588554},
    {18.770432, 18.770432, 17.332857, 18.770432, 17.332857, 16.593273, 17.332857, 16.593273, 17.332857, 16.593273, 17.332857, 18.770432, 17.332857, 18.770432, 18.770432},
    {18.936631, 18.936631, 17.619156, 18.936631, 17.619156, 16.954038, 17.619156, 16.954038, 17.619156, 16.954038, 17.619156, 18.936631, 17.619156, 18.936631, 18.936631},
    {19.08804, 19.08804, 17.876244, 19.08804, 17.876244, 17.274948, 17.876244, 17.274948, 17.876244, 17.274948, 17.876244, 19.08804, 17.876244, 19.08804, 19.08804},
    {19.225576, 19.225576, 18.107399, 19.225576, 18.107399, 17.561202, 18.107399, 17.561202, 18.107399, 17.561202, 18.107399, 19.225576, 18.107399, 19.225576, 19.225576},
    {19.350174, 19.350174, 18.315414, 19.350174, 18.315414, 17.817113, 18.315414, 17.817113, 18.315414, 17.817113, 18.315414, 19.350174, 18.315414, 19.350174, 19.350174},
    {19.462771, 19.462771, 18.502702, 19.462771, 18.502702, 18.046293, 18.502702, 18.046293, 18.502702, 18.046293, 18.502702, 19.462771, 18.502702, 19.462771, 19.462771}};
    // clang-format on

    EXPECT_EQ(temperatureSolution.size(), correctTemperatureSolution.size());

    for(auto i = 0u; i < temperatureSolution.size(); ++i)
    {
        for(auto j = 0u; j < temperatureSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctTemperatureSolution[i][j], temperatureSolution[i][j], 1e-6);
        }
    }
}
