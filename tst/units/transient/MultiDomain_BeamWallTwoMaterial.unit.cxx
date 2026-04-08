#include <cmath>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

#include "BeamBuilder.hxx"
#include "PrintSolution.hxx"
#include "TestMaterials.hxx"

//! Beam-shaped wall section with two materials: an exterior stucco layer
//! and a fiberglass-batt insulation layer. Initialised at high humidity
//! (RH = 0.9999) and warm temperature, with a colder/drier interior
//! boundary. Intended as a reproducer for solver divergence issues
//! reported on multi-material beam structures at near-saturation states.
TEST(MultiDomain_BeamWall_StuccoFiberglass, TwoMaterialHighHumidity)
{
    HygroThermFEM::MultiDomain multiDomain;

    // Initial state: warm, near-saturation (the regime where the solver
    // currently struggles).
    constexpr HygroThermFEM::State initialState({.temperature = 30.0,
                                                 .humidity = 0.9999,
                                                 .pressure = 101325.0,
                                                 .liquidPercent = 1.0});

    // Register both materials so they are available to the elements.
    const auto & stucco =
      multiDomain.materials().createSolidMaterial(TestHelper::Stucco());
    const auto & fiberglass =
      multiDomain.materials().createSolidMaterial(TestHelper::FiberglassBatts());

    // Beam geometry:
    //   - 0.02 m of exterior stucco (left), 4 elements wide
    //   - 0.10 m of fiberglass batt insulation (right), 10 elements wide
    //   - Total height 0.05 m, 4 element rows
    TestHelper::BeamBuilder builder(multiDomain);
    builder.xStart(0.0)
      .height(0.05)
      .numElementsY(4)
      .state(initialState)
      .addSegment({.material = stucco.name(), .numElementsX = 4, .width = 0.02})
      .addSegment({.material = fiberglass.name(), .numElementsX = 10, .width = 0.10})
      .build();

    // Boundary conditions:
    //   - Left edge (exterior stucco): cold, fully saturated air
    //   - Right edge (interior fiberglass): warmer, drier air
    constexpr auto hc = 10.0;
    const HygroThermFEM::FixedBCHCCoefficients exteriorBc{
      /*airTemperature=*/0.0, hc, /*airHumidity=*/1.0};
    const HygroThermFEM::FixedBCHCCoefficients interiorBc{
      /*airTemperature=*/20.0, hc, /*airHumidity=*/0.5};

    for(auto [i1, i2] : builder.leftEdge())
    {
        multiDomain.createBC_FixedHc(i1, i2, exteriorBc);
    }
    for(auto [i1, i2] : builder.rightEdge())
    {
        multiDomain.createBC_FixedHc(i1, i2, interiorBc);
    }

    // Drive a few transient timesteps and verify the solver does not
    // explode (no NaN/Inf, values stay within physical bounds).
    constexpr double dTime = 3600.0;
    constexpr int nSteps = 3;

    auto temperatures =
      multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
    auto humidities =
      multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);

    for(int step = 0; step < nSteps; ++step)
    {
        const auto solution = multiDomain.transient(
          temperatures, humidities, dTime, static_cast<size_t>(step));

        for(const double t : solution.temperature)
        {
            EXPECT_TRUE(std::isfinite(t)) << "Temperature blew up at step " << step;
            EXPECT_GT(t, -50.0) << "Temperature unphysically low at step " << step;
            EXPECT_LT(t, 200.0) << "Temperature unphysically high at step " << step;
        }
        for(const double h : solution.humidity)
        {
            EXPECT_TRUE(std::isfinite(h)) << "Humidity blew up at step " << step;
            EXPECT_GE(h, 0.0) << "Humidity below 0 at step " << step;
            EXPECT_LE(h, 1.0 + 1e-9) << "Humidity above 1 at step " << step;
        }

        temperatures = solution.temperature;
        humidities = solution.humidity;
    }
}
