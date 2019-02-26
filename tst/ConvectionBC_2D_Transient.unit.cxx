#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::NodePool;
using HygroThermFEM::MaterialPool;

class ConvectionBC_2D_Transient : public testing::Test
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

TEST_F(ConvectionBC_2D_Transient, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Three elements with simple convection BC.");

    NodePool::Instance().createNode(1, 0.2, 0.05);
    NodePool::Instance().createNode(2, 0.2, 0.00);
    NodePool::Instance().createNode(3, 0.1, 0.05);
    NodePool::Instance().createNode(4, 0.1, 0.00);
    NodePool::Instance().createNode(5, 0.0, 0.05);
    NodePool::Instance().createNode(6, 0.0, 0.00);

    auto & material = MaterialPool::Instance().createMaterial(
      "Test Material",
      2050,                       /// Density
      0.00,                       /// Porosity
      850,                        /// Specific Heat Capacity (dry)
      15E-6,                      /// Diffusion Resistance Factor
      {{0.0, 1.0}, {180, 1.0}},   /// Thermal Conductivity (as function of water content)
      {{0, 0},   /// Liquid Transportation Coefficient (as function of water content)
       {27, 1E-8},
       {45, 1.1E-8},
       {90, 2E-8},
       {126, 3.5E-8},
       {144, 5E-8},
       {162, 1E-7},
       {171, 2E-7},
       {180, 7E-7}},
      {{0, 0},   /// Moisture Storage Function (Water content as function of relative humidity)
       {0.5, 5.3},
       {0.65, 8.4},
       {0.8, 12},
       {0.93, 17},
       {0.95, 25},
       {0.99, 63},
       {0.995, 83},
       {0.999, 120},
       {1, 180}});

    HygroThermFEM::ThermalDomain domain;

    domain.createElement(3, 4, 2, 1, material.name());
    domain.createElement(6, 4, 3, 5, material.name());

    // Create Boundary Conditions
    const auto hc1 = 2.4;
    const auto temperatureAir1 = 20.0;

    const auto hc2 = 15.0;
    const auto temperatureAir2 = -18.0;

    domain.createConvectionBCFixedHc(1, 2, temperatureAir1, hc1);
    domain.createConvectionBCFixedHc(6, 5, temperatureAir2, hc2);

    const auto dTime = 3600;
    const auto nSteps = 4;

    auto temperatures = NodePool::Instance().properties(HygroThermFEM::Variable::temperature);
    std::vector<std::vector<double>> temperaturesSolution;
    std::vector<std::vector<HygroThermFEM::NodeFlux>> fluxSolution;

    for(unsigned i = 0; i < nSteps; ++i)
    {
        temperatures = domain.transient(temperatures, dTime).solution;
        temperaturesSolution.push_back(temperatures);
        fluxSolution.push_back(domain.flux());
    }

    std::vector<std::vector<double>> correctTemperatureSolution{
      {1.132130505, 1.132130505, -0.656245111, -0.656245111, -5.621029352, -5.621029352},
      {1.657776414, 1.657776414, -1.47222114, -1.47222114, -8.551769334, -8.551769334},
      {1.788809043, 1.788809043, -2.264759626, -2.264759626, -10.15443472, -10.15443472},
      {1.680631717, 1.680631717, -2.977820826, -2.977820826, -11.08768765, -11.08768765}};

    std::vector<std::vector<HygroThermFEM::NodeFlux>> correctFluxSolution{{{-17.88375616, 0},
                                                                           {-17.88375616, 0},
                                                                           {-33.76579929, 0},
                                                                           {-33.76579929, 0},
                                                                           {-49.64784241, 0},
                                                                           {-49.64784241, 0}},
                                                                          {{-31.29997554, 0},
                                                                           {-31.29997554, 0},
                                                                           {-51.04772874, 0},
                                                                           {-51.04772874, 0},
                                                                           {-70.79548195, 0},
                                                                           {-70.79548195, 0}},
                                                                          {{-40.53568669, 0},
                                                                           {-40.53568669, 0},
                                                                           {-59.71621882, 0},
                                                                           {-59.71621882, 0},
                                                                           {-78.89675094, 0},
                                                                           {-78.89675094, 0}},
                                                                          {{-46.58452543, 0},
                                                                           {-46.58452543, 0},
                                                                           {-63.84159682, 0},
                                                                           {-63.84159682, 0},
                                                                           {-81.09866822, 0},
                                                                           {-81.09866822, 0}}};

    EXPECT_EQ(temperaturesSolution.size(), correctTemperatureSolution.size());

    for(auto i = 0u; i < correctTemperatureSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctTemperatureSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctTemperatureSolution[i][j], temperaturesSolution[i][j], 1e-6);
            EXPECT_NEAR(correctFluxSolution[i][j].x, fluxSolution[i][j].x, 1e-6);
            EXPECT_NEAR(correctFluxSolution[i][j].y, fluxSolution[i][j].y, 1e-6);
        }
    }
}
