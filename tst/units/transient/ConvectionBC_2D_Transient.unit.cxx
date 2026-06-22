#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "TestMaterials.hxx"

TEST(ConvectionBC_2D_Transient, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Three elements with simple convection BC.");

    HygroThermFEM::MultiDomain multiDomain({.performThermal = true, .performMoisture = false});

    multiDomain.nodes().createNode({.index = 1, .x = 0.2, .y = 0.05});
    multiDomain.nodes().createNode({.index = 2, .x = 0.2, .y = 0.00});
    multiDomain.nodes().createNode({.index = 3, .x = 0.1, .y = 0.05});
    multiDomain.nodes().createNode({.index = 4, .x = 0.1, .y = 0.00});
    multiDomain.nodes().createNode({.index = 5, .x = 0.0, .y = 0.05});
    multiDomain.nodes().createNode({.index = 6, .x = 0.0, .y = 0.00});

    const auto & material = multiDomain.materials().createSolidMaterial(TestHelper::TestMaterial());

    multiDomain.createElement({.node1 = 3, .node2 = 4, .node3 = 2, .node4 = 1, .material = material.name()});
    multiDomain.createElement({.node1 = 6, .node2 = 4, .node3 = 3, .node4 = 5, .material = material.name()});

    // Create Boundary Conditions
    constexpr auto hc1 = 2.4;
    constexpr auto temperatureAir1 = 20.0;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff1{temperatureAir1, hc1};

    constexpr auto hc2 = 15.0;
    constexpr auto temperatureAir2 = -18.0;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff2{temperatureAir2, hc2};

    multiDomain.thermal().createBC_FixedHc(1, 2, bcCoeff1);
    multiDomain.thermal().createBC_FixedHc(6, 5, bcCoeff2);

    constexpr auto dTime = 3600;
    constexpr auto nSteps = 4;

    auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
    std::vector<std::vector<double>> temperaturesSolution;
    std::vector<std::vector<HygroThermFEM::NodeFlux>> fluxSolution;

    for(unsigned i = 0; i < nSteps; ++i)
    {
        temperatures = multiDomain.thermal().transient(temperatures, dTime).value().solution;
        temperaturesSolution.push_back(temperatures);
        fluxSolution.push_back(multiDomain.thermal().flux());
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
