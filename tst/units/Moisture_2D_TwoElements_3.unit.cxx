#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "TestMaterials.hxx"

using HygroThermFEM::State;

TEST(Moisture_2D_TwoElements_3, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Simple two elements example with moisture transfer.");

    HygroThermFEM::MultiDomain multiDomain({.performThermal = false, .performMoisture = true});

    const State state({
        .temperature = 20.0,
        .humidity = 0.0,
        .pressure = 101325.0,
        .liquidPercent = 1.0
    });

    multiDomain.nodes().createNode({.index = 1, .x = 0.15, .y = 0.05, .state = state});
    multiDomain.nodes().createNode({.index = 2, .x = 0.15, .y = 0, .state = state});
    multiDomain.nodes().createNode({.index = 3, .x = 0.05, .y = 0.05, .state = state});
    multiDomain.nodes().createNode({.index = 4, .x = 0.05, .y = 0, .state = state});
    multiDomain.nodes().createNode({.index = 5, .x = 0, .y = 0.05, .state = state});
    multiDomain.nodes().createNode({.index = 6, .x = 0, .y = 0, .state = state});

    const auto & material = multiDomain.materials().createSolidMaterial(TestHelper::CottaerSandstone());

    /// Create elements
    multiDomain.createElement({.node1 = 3, .node2 = 4, .node3 = 2, .node4 = 1, .material = material.name()});
    multiDomain.createElement({.node1 = 6, .node2 = 4, .node3 = 3, .node4 = 5, .material = material.name()});

    // Create Boundary Conditions
    constexpr auto airTemperature = 20.0;
    constexpr auto airHumidity = 0.5;

    const HygroThermFEM::TemperatureAndHumidity bcCoeff{airTemperature, airHumidity};

    multiDomain.moisture().createBC_FixedHumidity(5, 6, bcCoeff);

    constexpr auto dTime = 3600;
    constexpr auto nSteps = 24;

    auto humidities = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);
    std::vector<double> timesteps;
    std::vector<std::vector<double>> humiditySolution;
    std::vector<std::vector<HygroThermFEM::NodeFlux>> fluxSolution;

    for(unsigned i = 0; i < nSteps; ++i)
    {
        auto solution = multiDomain.moisture().transient(humidities, dTime);
        humidities = solution.solution;
        timesteps.push_back(solution.dTime);
        auto humidityContent = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);
        humiditySolution.push_back(humidityContent);
        fluxSolution.push_back(multiDomain.moisture().flux());
    }

    std::vector<double> correctTimesteps{3600, 3600, 3600, 3600, 3600, 3600, 3600, 3600,
                                         3600, 3600, 3600, 3600, 3600, 3600, 3600, 3600,
                                         3600, 3600, 3600, 3600, 3600, 3600, 3600, 3600};

    std::vector<std::vector<double>> correctHumiditySolution{
      {2e-06, 2e-06, 0.001245, 0.001245, 0.5, 0.5},
      {7e-06, 7e-06, 0.002486, 0.002486, 0.5, 0.5},
      {1.4e-05, 1.4e-05, 0.003722, 0.003722, 0.5, 0.5},
      {2.4e-05, 2.4e-05, 0.004953, 0.004953, 0.5, 0.5},
      {3.6e-05, 3.6e-05, 0.00618, 0.00618, 0.5, 0.5},
      {5e-05, 5e-05, 0.007402, 0.007402, 0.5, 0.5},
      {6.7e-05, 6.7e-05, 0.00862, 0.00862, 0.5, 0.5},
      {8.6e-05, 8.6e-05, 0.009833, 0.009833, 0.5, 0.5},
      {0.000107, 0.000107, 0.011042, 0.011042, 0.5, 0.5},
      {0.000131, 0.000131, 0.012246, 0.012246, 0.5, 0.5},
      {0.000157, 0.000157, 0.013446, 0.013446, 0.5, 0.5},
      {0.000185, 0.000185, 0.014642, 0.014642, 0.5, 0.5},
      {0.000215, 0.000215, 0.015833, 0.015833, 0.5, 0.5},
      {0.000248, 0.000248, 0.017019, 0.017019, 0.5, 0.5},
      {0.000283, 0.000283, 0.018201, 0.018201, 0.5, 0.5},
      {0.00032, 0.00032, 0.019379, 0.019379, 0.5, 0.5},
      {0.00036, 0.00036, 0.020552, 0.020552, 0.5, 0.5},
      {0.000401, 0.000401, 0.021721, 0.021721, 0.5, 0.5},
      {0.000445, 0.000445, 0.022886, 0.022886, 0.5, 0.5},
      {0.000491, 0.000491, 0.024047, 0.024047, 0.5, 0.5},
      {0.000539, 0.000539, 0.025203, 0.025203, 0.5, 0.5},
      {0.000589, 0.000589, 0.026355, 0.026355, 0.5, 0.5},
      {0.000642, 0.000642, 0.027503, 0.027503, 0.5, 0.5},
      {0.000696, 0.000696, 0.028646, 0.028646, 0.5, 0.5}};


    EXPECT_EQ(correctTimesteps.size(), timesteps.size());

    for(size_t i = 0u; i < correctTimesteps.size(); ++i)
    {
        EXPECT_NEAR(correctTimesteps[i], timesteps[i], 1e-6);
    }

    EXPECT_EQ(correctHumiditySolution.size(), humiditySolution.size());

    for(auto i = 0u; i < humiditySolution.size(); ++i)
    {
        for(auto j = 0u; j < humiditySolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctHumiditySolution[i][j], humiditySolution[i][j], 1e-6);
        }
    }
}
