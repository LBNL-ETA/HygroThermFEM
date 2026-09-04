#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

#include "BeamBuilder.hxx"
#include "SteadyVsTransientHarness.hxx"
#include "TestMaterials.hxx"

// Steady state against transient on a two-material wall, once for each moisture modelling
// option. What these tests are for is described in SteadyVsTransientHarness.hxx; this file
// supplies the wall, and each test states the options it ticks, how long it marches and how
// far the two solvers are allowed to sit apart.
//
// The wall is the same assembly as ThermSample_StuccoFiberglass: Stucco 76.3 mm on the
// exterior side, Fiberglass Batts 80 mm on the interior side, as a beam of 8 + 8 elements,
// one element tall. Constant films: -18 C / RH 0.5 / hc 26 W/m2K on the stucco face,
// 21 C / RH 0.5 / hc 4.65 W/m2K on the fiberglass face, initial state 21 C / RH 0.5
// throughout. Every case uses that same wall and those same films, so the only difference
// between cases is the option under test.
//
// With no vapour barrier in it this wall condenses at the stucco/fiberglass interface, and
// its thick stucco layer holds enough moisture that ten simulated years still leave the
// transient short of rest. That, not solver error, is what sets the humidity tolerances
// below; lengthening the march is the lever if they ever need to be tighter. The
// four-material wall in SteadyVsTransient_StuccoWall settles far closer and is the sharper
// check.
namespace
{
    using namespace TestHelper::SteadyVsTransient;

    void buildBeam(HygroThermFEM::MultiDomain & multiDomain)
    {
        const auto & stucco = multiDomain.materials().createSolidMaterial(TestHelper::Stucco());
        const auto & fiberglass =
          multiDomain.materials().createSolidMaterial(TestHelper::FiberglassBatts());

        TestHelper::BeamBuilder builder(multiDomain);
        builder.xStart(0.0)
          .height(0.0728)
          .numElementsY(1)
          .state({.temperature = 21.0, .humidity = 0.5, .pressure = 101325.0, .liquidPercent = 1.0})
          .addSegment({.material = stucco.name(), .numElementsX = 8, .width = 0.0763})
          .addSegment({.material = fiberglass.name(), .numElementsX = 8, .width = 0.08})
          .build();

        const HygroThermFEM::FixedBCHCCoefficients exterior{-18.0, 26.0, 0.5};
        const HygroThermFEM::FixedBCHCCoefficients interior{21.0, 4.65, 0.5};
        builder.applyBC_FixedHc(TestHelper::BeamBuilder::Edge::Left, exterior);
        builder.applyBC_FixedHc(TestHelper::BeamBuilder::Edge::Right, interior);
    }
}   // namespace

TEST(SteadyVsTransient_StuccoFiberglassBeam, AllOff)
{
    expectSteadyMatchesTransient(
      buildBeam,
      {.name = "all options off",
       // Measured 2.06e-4 in humidity at the saturation front; temperature exact.
       .tolerances = {.temperature = 1e-12, .humidity = 2.1e-4}});
}

TEST(SteadyVsTransient_StuccoFiberglassBeam, LiquidTransportOnly)
{
    expectSteadyMatchesTransient(
      buildBeam,
      {.name = "water liquid transportation only",
       .options = {.liquidTransport = true},
       // Measured 7.91e-12 in humidity; temperature exact. Liquid redistribution brings this
       // wall properly to rest, so it is by far the best behaved case here.
       .tolerances = {.temperature = 1e-12, .humidity = 8e-12}});
}

TEST(SteadyVsTransient_StuccoFiberglassBeam, HeatOfEvaporationOnly)
{
    expectSteadyMatchesTransient(
      buildBeam,
      {.name = "heat of evaporation only",
       .options = {.heatOfEvaporation = true},
       // Measured 9.04e-7 in temperature and 1.64e-4 in humidity. The only case here whose
       // temperature is not exact, being the only one that couples moisture back into heat.
       .tolerances = {.temperature = 9.1e-7, .humidity = 1.65e-4}});
}

//! Capillary conduction is heat carried by the liquid flux, and with liquid transport off
//! there is no such flux, so the term is inactive here and this case behaves as "all off".
//! (Until 2026-09-03 it ran on a fictitious flux built from the vapour-only humidity
//! gradients, and the transient blew up to -273 C mid-wall at 5 h steps.)
TEST(SteadyVsTransient_StuccoFiberglassBeam, CapillaryConductionOnly)
{
    expectSteadyMatchesTransient(
      buildBeam,
      {.name = "capillary conduction only",
       .options = {.capillaryConduction = true},
       // Measured 2.06e-4 in humidity, identical to "all options off" to every digit, which
       // is the term being inert without a liquid flux.
       .tolerances = {.temperature = 1e-12, .humidity = 2.1e-4}});
}

//! Liquid transport plus capillary conduction: the pair that actually exercises the term.
TEST(SteadyVsTransient_StuccoFiberglassBeam, LiquidTransportAndCapillaryConduction)
{
    expectSteadyMatchesTransient(
      buildBeam,
      {.name = "water liquid transportation + capillary conduction",
       .options = {.liquidTransport = true, .capillaryConduction = true},
       // Measured 7.83e-12 in humidity; temperature exact, as with liquid transport alone.
       .tolerances = {.temperature = 1e-12, .humidity = 8e-12}});
}

TEST(SteadyVsTransient_StuccoFiberglassBeam, VaporDiffusionOnly)
{
    expectSteadyMatchesTransient(
      buildBeam,
      {.name = "vapor diffusion only",
       .options = {.vaporDiffusionConduction = true},
       // Measured 5.84e-10 in temperature and 2.06e-4 in humidity.
       .tolerances = {.temperature = 6e-10, .humidity = 2.1e-4}});
}
