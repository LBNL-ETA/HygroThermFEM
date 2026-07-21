#include <array>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "SlabCreator.hxx"
#include "TestMaterials.hxx"

/////////////////////////////////////////////////////////////////////////////////////
/// TOPAZ2D Validation Status Report (Drayton, UCRL-ID-106550, December 1990),
/// section 3.9 -- Plane Transient Linear, Spatial Element Generation.
///
/// Unit slab (alpha = 1, k = 1) at zero initial temperature, adiabatic at x = 0,
/// surface at x = 1 held at 0, spatially varying generation q(x) = 1 - (x/L)^2
/// applied per element at its midpoint. Analytical solution: the odd-cosine
/// eigenfunction expansion (hygrothermfem_python, analytic.slab_generation_spatial;
/// its docstring records why the report's printed series, which also carries
/// even-index terms, does not satisfy T(L) = 0). Expected values are the exact
/// odd-index series at the report's checkpoints; the tolerance is the measured
/// backward-Euler discretization band at dt = 0.01 on ten elements (max 2.5e-3).
/////////////////////////////////////////////////////////////////////////////////////

TEST(Topaz2DReport_GenerationSpatial, ParabolicGeneration)
{
    SCOPED_TRACE("Begin Test: unit slab, spatially varying generation.");

    HygroThermFEM::MultiDomain multiDomain({.performThermal = true, .performMoisture = false});

    constexpr HygroThermFEM::State state({
        .temperature = 0.0,
        .humidity = 0.0,
        .pressure = 101325.0,
        .liquidPercent = 0
    });

    auto params = TestHelper::TestMaterial();
    params.density = 1.0;
    params.heatCapacity = 1.0;
    const auto & material = multiDomain.materials().createSolidMaterial(params);

    TestHelper::SlabBuilder(multiDomain)
        .gridXCoordinates({0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0})
        .height(0.05)
        .material(material.name())
        .state(state)
        .build();

    constexpr auto tSurface = 0.0;
    multiDomain.thermal().createBC_FixedTemperature(21, 22, tSurface);

    // q(x) = 1 - x^2 evaluated at element midpoints, in element creation order
    // (SlabBuilder builds columns left to right).
    std::vector<double> generation;
    generation.reserve(10);
    for(std::size_t elem = 0; elem < 10; ++elem)
    {
        const double midpoint = 0.1 * static_cast<double>(elem) + 0.05;
        generation.push_back(1.0 - midpoint * midpoint);
    }
    multiDomain.thermal().setVolumetricSource(generation);

    constexpr auto dTime = 0.01;
    constexpr auto nSteps = 200;

    const auto solution = multiDomain.thermal().transientMultiStep(
      HygroThermFEM::Variable::temperature, dTime, nSteps);

    for(std::size_t step = 0; step < solution.size(); ++step)
    {
    }

    // Exact odd-index series at x = 0, 0.5, 1.0 for steps 10, 50, 100, 200.
    const std::vector<std::pair<std::size_t, std::array<double, 3>>> checkpoints{
      {10, {0.090037642, 0.065649588, 0.0}},
      {50, {0.294859968, 0.210744621, 0.0}},
      {100, {0.381194973, 0.271792725, 0.0}},
      {200, {0.413658491, 0.294747898, 0.0}},
    };

    for(const auto & [step, expected] : checkpoints)
    {
        for(std::size_t pos = 0; pos < expected.size(); ++pos)
        {
            EXPECT_NEAR(expected[pos], solution[step - 1][pos * 10], 0.005)
              << "step " << step << ", position " << pos;
        }
    }
}
