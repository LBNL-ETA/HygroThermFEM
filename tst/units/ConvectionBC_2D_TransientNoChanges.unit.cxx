#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

TEST(ConvectionBC_2D_TransientNoChanges, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Three elements with simple convection BC.");

    HygroThermFEM::MultiDomain multiDomain({.performThermal = true, .performMoisture = false});

    std::vector<double> gridXCoordinates{0, 0.05, 0.1, 0.15};

    const HygroThermFEM::State state({
        .temperature = 20.0,
        .humidity = 0.0,
        .pressure = 101325.0,
        .liquidPercent = 0
    });
    for(auto val : gridXCoordinates)
    {
        multiDomain.nodes().createNode({.x = val, .y = 0.00, .state = state});
        multiDomain.nodes().createNode({.x = val, .y = 0.05, .state = state});
    }

    // Material Properties (Cottaer Sandstone, using C++20 designated initializers)
    const auto & material = multiDomain.materials().createSolidMaterial({
        .name = "Cottaer Sandstone",
        .thermalConductivityDry = 1.8,
        .density = 2050.0,
        .porosity = 0.22,
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

    /// Create elements
    for(size_t idx = 1u; idx <= (multiDomain.nodes().maxIndex() - 2) / 2; ++idx)
    {
        const auto index1 = 2u * idx + 1u;
        const auto index2 = 2u * idx + 2u;
        const auto index3 = 2u * idx;
        const auto index4 = 2u * idx - 1u;
        multiDomain.createElement({.node1 = index1, .node2 = index2, .node3 = index3, .node4 = index4, .material = material.name()});
    }

    // Create Boundary Conditions
    constexpr auto tSurface = 20.0;
    constexpr auto hc = 1.0;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff{tSurface, hc};

    multiDomain.thermal().createBC_FixedHc(1, 2, bcCoeff);

    constexpr auto dTime = 36000;
    constexpr auto nSteps = 4;

    auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
    std::vector<std::vector<double>> solution;

    for(unsigned i = 0; i < nSteps; ++i)
    {
        temperatures = multiDomain.thermal().transient(temperatures, dTime).solution;
        solution.push_back(temperatures);
    }

    std::vector<std::vector<double>> correctSolution{
      {20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0},
      {20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0},
      {20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0},
      {20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0}};

    EXPECT_EQ(solution.size(), correctSolution.size());

    for(auto i = 0u; i < correctSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctSolution[i][j], solution[i][j], 1e-6);
        }
    }
}
