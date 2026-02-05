#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

/////////////////////////////////////////////////////////////////////////////////////
/// Transient heat transfer example on Sandstone specimen using data from database
///   Lumped mass matrix
///   Time-step 1 hour
///   Six nodes block at initial temperatures in nodes of 0 degrees Celsius
///   Initial Black Body Radiation boundary condition at nodes 5 and 6
///   Solution achieved with nonlinear solver
/////////////////////////////////////////////////////////////////////////////////////

TEST(BlackBodyBC_2D_1, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Two elementsCreator example with transient.");

    HygroThermFEM::MultiDomain multiDomain({.performThermal = true, .performMoisture = false});

    // Enter nodes. Arguments are: node number, x-coordinate, y-coordinate, initial temperature

    // same temperature in every node (humidity and pressure irrelevant for this example)
    const auto state = HygroThermFEM::State{
        .temperature = 0,
        .humidity = 0,
        .pressure = 101325,
        .liquidPercent = 0
    };

    multiDomain.nodes().createNode({.index = 1, .x = 0.15, .y = 0.05, .state = state});
    multiDomain.nodes().createNode({.index = 2, .x = 0.15, .y = 0, .state = state});
    multiDomain.nodes().createNode({.index = 3, .x = 0.05, .y = 0.05, .state = state});
    multiDomain.nodes().createNode({.index = 4, .x = 0.05, .y = 0, .state = state});
    multiDomain.nodes().createNode({.index = 5, .x = 0, .y = 0.05, .state = state});
    multiDomain.nodes().createNode({.index = 6, .x = 0, .y = 0, .state = state});

    // Material Properties (Cottaer Sandstone - non porous, using C++20 designated initializers)
    const auto & material = multiDomain.materials().createSolidMaterial({
        .name = "Cottaer Sandstone - non porous",
        .thermalConductivityDry = 1.8,
        .density = 2050.0,
        .porosity = 0.0,
        .heatCapacity = 850.0,
        .diffusionResistanceFactor = 15.0,
        .thermalConductivityMoistureDependent = {{0.0, 1.8}, {180, 1.8}},
        .moistureDependentMeasurementTemperature = 0.0,
        .thermalConductivityTemperatureDependent = {{0.0, 1.8}, {1, 1.8}},
        .temperatureDependentMeasurementHumidity = 0.0,
        .liquidTransportCurve = {{0, 0}, {27, 1E-8}, {45, 1.1E-8}, {90, 2E-8}, {126, 3.5E-8},
                                 {144, 5E-8}, {162, 1E-7}, {171, 2E-7}, {180, 7E-7}},
        .sorptionCurve = {{0, 0}, {0.5, 5.3}, {0.65, 8.4}, {0.8, 12}, {0.93, 17},
                          {0.95, 25}, {0.99, 63}, {0.995, 83}, {0.999, 120}, {1, 180}}
    });

    multiDomain.createElement({.node1 = 3, .node2 = 4, .node3 = 2, .node4 = 1, .material = material.name()});
    multiDomain.createElement({.node1 = 6, .node2 = 4, .node3 = 3, .node4 = 5, .material = material.name()});

    // Create Boundary Conditions
    constexpr auto tRadiation = 20.0;
    constexpr auto surfaceEmissivity = 0.84;

    multiDomain.thermal().createBC_BlackBodyRadiation(5, 6, surfaceEmissivity, tRadiation);

    constexpr auto dTime = 3600;
    constexpr auto nSteps = 4;

    auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
    std::vector<std::vector<double>> solution;

    for(unsigned i = 0; i < nSteps; ++i)
    {
        temperatures = multiDomain.thermal().transient(temperatures, dTime).solution;
        solution.push_back(temperatures);
    }

    std::vector<std::vector<double>> correctSolution = {
      {0.442761, 0.442761, 1.038062, 1.038062, 2.382485, 2.382485},
      {1.155794, 1.155794, 2.114481, 2.114482, 3.679278, 3.679278},
      {1.987435, 1.987435, 3.105594, 3.105594, 4.664103, 4.664103},
      {2.853862, 2.853862, 4.018792, 4.018792, 5.522119, 5.522119}};

    EXPECT_EQ(solution.size(), correctSolution.size());

    for(auto i = 0u; i < solution.size(); ++i)
    {
        for(auto j = 0u; j < solution[i].size(); ++j)
        {
            EXPECT_NEAR(correctSolution[i][j], solution[i][j], 1e-6);
        }
    }
}
