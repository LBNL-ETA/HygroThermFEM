#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

#include "BeamBuilder.hxx"
#include "SteadyVsTransientHarness.hxx"
#include "TestMaterials.hxx"

// Steady state against transient on a four-material wall containing a vapour barrier. What
// these tests are for is described in SteadyVsTransientHarness.hxx; this file supplies the
// wall, and each test states the options it ticks, how long it marches and how far the two
// solvers are allowed to sit apart.
//
// The wall is the THERM sample "Stucco Wall - Moisture" as a beam, with THERM's own layer
// split: Stucco 25.4 mm / 3 elements | Laminated panel 19.05 mm / 2 | Fiberglass Batts
// 76.2 mm / 7 | Gypsum Board Interior 25.4 mm / 3, one element tall. It carries the same
// constant films as the two-material wall so the two are directly comparable: -18 C / RH 0.5
// / hc 26 W/m2K on the stucco face, 21 C / RH 0.5 / hc 4.65 W/m2K on the gypsum face,
// initial state 21 C / RH 0.5 throughout.
//
// WHY THIS WALL IS HERE and not just more of the same. The laminated panel is a near vapour
// barrier, with a diffusion resistance factor two orders above the fiberglass beside it, and
// that changes both the physics and the numbers:
//
//  - The condensation plane pins BEHIND the barrier, at the stucco/panel face, rather than
//    at the single interface of the two-material wall. Almost the whole vapour pressure drop
//    sits across the panel.
//  - Switching liquid transport on wets the stucco from about 0.50 to 0.94 while the
//    interior of the barrier does not move at all: liquid bypasses a barrier that vapour
//    cannot cross. The two-material wall has no barrier and cannot exercise this.
//  - Its stucco layer is 25 mm rather than 76 mm, so it holds far less moisture and settles
//    much closer to rest in the same ten simulated years. That is what lets the tolerances
//    here be orders tighter than the two-material wall can support.
//
// The option matrix itself is covered by the two-material wall, so only the cases the
// barrier changes are repeated here: nothing on, liquid transport, heat of evaporation (the
// latent release now happens behind the barrier), and liquid with capillary conduction.
// Capillary alone is inert without a liquid flux, and vapour diffusion alone is a small
// perturbation of nothing-on; both are already established there.
namespace
{
    using namespace TestHelper::SteadyVsTransient;

    void buildWall(HygroThermFEM::MultiDomain & multiDomain)
    {
        const auto & stucco = multiDomain.materials().createSolidMaterial(TestHelper::Stucco());
        const auto & panel =
          multiDomain.materials().createSolidMaterial(TestHelper::LaminatedPanel());
        const auto & fiberglass =
          multiDomain.materials().createSolidMaterial(TestHelper::FiberglassBatts());
        const auto & gypsum =
          multiDomain.materials().createSolidMaterial(TestHelper::GypsumBoardInterior());

        TestHelper::BeamBuilder builder(multiDomain);
        builder.xStart(0.0)
          .height(0.0728)
          .numElementsY(1)
          .state({.temperature = 21.0, .humidity = 0.5, .pressure = 101325.0, .liquidPercent = 1.0})
          .addSegment({.material = stucco.name(), .numElementsX = 3, .width = 0.0254})
          .addSegment({.material = panel.name(), .numElementsX = 2, .width = 0.01905})
          .addSegment({.material = fiberglass.name(), .numElementsX = 7, .width = 0.0762})
          .addSegment({.material = gypsum.name(), .numElementsX = 3, .width = 0.0254})
          .build();

        const HygroThermFEM::FixedBCHCCoefficients exterior{-18.0, 26.0, 0.5};
        const HygroThermFEM::FixedBCHCCoefficients interior{21.0, 4.65, 0.5};
        builder.applyBC_FixedHc(TestHelper::BeamBuilder::Edge::Left, exterior);
        builder.applyBC_FixedHc(TestHelper::BeamBuilder::Edge::Right, interior);
    }
}   // namespace

TEST(SteadyVsTransient_StuccoWall, AllOff)
{
    expectSteadyMatchesTransient(
      buildWall,
      {.name = "all options off",
       // Measured 9.79e-6 in humidity at the barrier face; temperature exact. Twenty times
       // closer than the same case on the two-material wall.
       .tolerances = {.temperature = 1e-11, .humidity = 1.2e-5}});
}

//! Liquid transport wets the stucco from 0.50 to 0.94 and leaves everything inside the
//! barrier untouched: the bypass this wall exists to cover.
TEST(SteadyVsTransient_StuccoWall, LiquidTransportOnly)
{
    expectSteadyMatchesTransient(
      buildWall,
      {.name = "water liquid transportation only",
       .options = {.liquidTransport = true},
       // Measured 8.79e-11 in humidity; temperature exact.
       .tolerances = {.temperature = 1e-11, .humidity = 2e-10}});
}

//! The latent release now happens behind the barrier, and the wall runs about 1 C warmer
//! through the fiberglass than with the option off.
TEST(SteadyVsTransient_StuccoWall, HeatOfEvaporationOnly)
{
    expectSteadyMatchesTransient(
      buildWall,
      {.name = "heat of evaporation only",
       .options = {.heatOfEvaporation = true},
       // Measured 2.85e-7 in temperature and 1.02e-4 in humidity: the loosest case on this
       // wall, and still tighter than the two-material wall manages with nothing on.
       .tolerances = {.temperature = 6e-7, .humidity = 1.2e-4}});
}

TEST(SteadyVsTransient_StuccoWall, LiquidTransportAndCapillaryConduction)
{
    expectSteadyMatchesTransient(
      buildWall,
      {.name = "water liquid transportation + capillary conduction",
       .options = {.liquidTransport = true, .capillaryConduction = true},
       // Measured 4.28e-12 in temperature and 8.73e-11 in humidity.
       .tolerances = {.temperature = 1e-11, .humidity = 2e-10}});
}
