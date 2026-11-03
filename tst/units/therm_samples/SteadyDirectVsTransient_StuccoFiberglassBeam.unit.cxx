#include <cmath>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

#include "BeamBuilder.hxx"
#include "SteadyVsTransientHarness.hxx"
#include "TestMaterials.hxx"

// The direct steady solve (no time variable) against the transient brought to rest, on the
// same two-material wall and the same six modelling cases as
// SteadyVsTransient_StuccoFiberglassBeam, which runs the pseudo-time march. The wall, the
// films and the schedule are described there. The tolerances below are the march's; the direct
// solve must earn the same ones, since both solve the same equations.
//
// The last test puts the two steady solvers side by side on the case that broke the previous
// direct solve, liquid transport near saturation.
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

TEST(SteadyDirectVsTransient_StuccoFiberglassBeam, AllOff)
{
    expectSteadyMatchesTransient(buildBeam,
                                 {.name = "all options off, direct",
                                  .tolerances = {.temperature = 1e-11, .humidity = 2.4e-4},
                                  .method = HygroThermFEM::SteadyMethod::direct});
}

TEST(SteadyDirectVsTransient_StuccoFiberglassBeam, LiquidTransportOnly)
{
    expectSteadyMatchesTransient(buildBeam,
                                 {.name = "water liquid transportation only, direct",
                                  .options = {.liquidTransport = true},
                                  .tolerances = {.temperature = 1e-11, .humidity = 2e-11},
                                  .method = HygroThermFEM::SteadyMethod::direct});
}

TEST(SteadyDirectVsTransient_StuccoFiberglassBeam, HeatOfEvaporationOnly)
{
    expectSteadyMatchesTransient(buildBeam,
                                 {.name = "heat of evaporation only, direct",
                                  .options = {.heatOfEvaporation = true},
                                  .tolerances = {.temperature = 2e-6, .humidity = 1.9e-4},
                                  .method = HygroThermFEM::SteadyMethod::direct});
}

TEST(SteadyDirectVsTransient_StuccoFiberglassBeam, CapillaryConductionOnly)
{
    expectSteadyMatchesTransient(buildBeam,
                                 {.name = "capillary conduction only, direct",
                                  .options = {.capillaryConduction = true},
                                  .tolerances = {.temperature = 1e-11, .humidity = 2.4e-4},
                                  .method = HygroThermFEM::SteadyMethod::direct});
}

TEST(SteadyDirectVsTransient_StuccoFiberglassBeam, LiquidTransportAndCapillaryConduction)
{
    expectSteadyMatchesTransient(
      buildBeam,
      {.name = "water liquid transportation + capillary conduction, direct",
       .options = {.liquidTransport = true, .capillaryConduction = true},
       .tolerances = {.temperature = 1e-11, .humidity = 2e-11},
       .method = HygroThermFEM::SteadyMethod::direct});
}

TEST(SteadyDirectVsTransient_StuccoFiberglassBeam, VaporDiffusionOnly)
{
    expectSteadyMatchesTransient(buildBeam,
                                 {.name = "vapor diffusion only, direct",
                                  .options = {.vaporDiffusionConduction = true},
                                  .tolerances = {.temperature = 1.5e-9, .humidity = 2.4e-4},
                                  .method = HygroThermFEM::SteadyMethod::direct});
}

//! The two steady solvers on the same wall with liquid transport on: the case on which the
//! previous direct solve oscillated. Both solve the same equations, so they must land on the
//! same field to solver precision, independently of the transient.
TEST(SteadyDirectVsTransient_StuccoFiberglassBeam, DirectAgreesWithMarch)
{
    const ModellingCase options{.name = "liquid transport, direct against march",
                                .options = {.liquidTransport = true}};
    ModellingCase direct{options};
    direct.method = HygroThermFEM::SteadyMethod::direct;

    const auto march = runSteady(buildBeam, options);
    const auto directRun = runSteady(buildBeam, direct);

    EXPECT_TRUE(march.converged) << "march stopped on its pass budget after " << march.passes;
    EXPECT_TRUE(directRun.converged)
      << "direct solve stopped on its pass budget after " << directRun.passes;

    constexpr double temperatureTolerance{1e-9};
    constexpr double humidityTolerance{1e-9};
    const auto temperatureDrift =
      worstDrift(directRun.fields.temperature, march.fields.temperature);
    const auto humidityDrift = worstDrift(directRun.fields.humidity, march.fields.humidity);
    EXPECT_LE(temperatureDrift.worst, temperatureTolerance)
      << "temperature differs by " << temperatureDrift.worst << " at column "
      << temperatureDrift.column + 1u;
    EXPECT_LE(humidityDrift.worst, humidityTolerance)
      << "humidity differs by " << humidityDrift.worst << " at column "
      << humidityDrift.column + 1u;
}
