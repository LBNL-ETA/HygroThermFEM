#include <stdexcept>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

#include "BeamBuilder.hxx"
#include "TestMaterials.hxx"

// The direct steady solve (MultiDomain::steadyState with SteadyMethod::direct) on the small
// cases the march is pinned on, plus the cases it must refuse.
namespace
{
    //! The two-element conduction strip of ConvectionBC_2D_SteadyState, with the same films.
    void buildConductionStrip(HygroThermFEM::MultiDomain & multiDomain)
    {
        multiDomain.nodes().createNode({.x = 15, .y = 5});
        multiDomain.nodes().createNode({.x = 15, .y = 0});
        multiDomain.nodes().createNode({.x = 5, .y = 5});
        multiDomain.nodes().createNode({.x = 5, .y = 0});
        multiDomain.nodes().createNode({.x = 0, .y = 5});
        multiDomain.nodes().createNode({.x = 0, .y = 0});

        const auto & material =
          multiDomain.materials().createSolidMaterial(TestHelper::TestMaterial());

        multiDomain.createElement(
          {.node1 = 3, .node2 = 4, .node3 = 2, .node4 = 1, .material = material.name()});
        multiDomain.createElement(
          {.node1 = 6, .node2 = 4, .node3 = 3, .node4 = 5, .material = material.name()});

        const HygroThermFEM::FixedBCHCCoefficients exterior{-18.0, 20.0, 0.0};
        const HygroThermFEM::FixedBCHCCoefficients interior{21.0, 2.4, 0.0};
        multiDomain.createBC_FixedHc(1, 2, exterior);
        multiDomain.createBC_FixedHc(6, 5, interior);
    }

    //! A sandstone strip with nothing applied to its faces by the caller.
    TestHelper::BeamBuilder buildSandstoneStrip(HygroThermFEM::MultiDomain & multiDomain)
    {
        const auto & sandstone =
          multiDomain.materials().createSolidMaterial(TestHelper::CottaerSandstone());
        TestHelper::BeamBuilder builder(multiDomain);
        builder.xStart(0.0)
          .height(0.01)
          .numElementsY(1)
          .state({.temperature = 20.0, .humidity = 0.5, .pressure = 101325.0, .liquidPercent = 1.0})
          .addSegment({.material = sandstone.name(), .numElementsX = 4, .width = 0.1})
          .build();
        return builder;
    }
}   // namespace

//! A linear thermal problem: one pass, one linear solve, the recorded field to solver precision.
TEST(SteadyStateDirect, LinearConductionMatchesRecordedField)
{
    HygroThermFEM::MultiDomain multiDomain({.performThermal = true, .performMoisture = false});
    buildConductionStrip(multiDomain);

    const auto solution = multiDomain.steadyState(HygroThermFEM::SteadyMethod::direct);

    EXPECT_TRUE(solution.converged);
    const std::vector<double> recorded{
      -17.87392241, -17.87392241, 7.341594828, 7.341594828, 19.94935345, 19.94935345};
    ASSERT_EQ(solution.temperature.size(), recorded.size());
    for(std::size_t node = 0u; node < recorded.size(); ++node)
    {
        EXPECT_NEAR(solution.temperature[node], recorded[node], 1e-6) << "node " << node + 1u;
    }
}

//! A moisture domain with no boundary exchange has no direct steady solution: the equations
//! fix the humidity profile only up to the total water. The solve says so instead of
//! factorising a singular system.
TEST(SteadyStateDirect, SealedMoistureDomainIsRefused)
{
    HygroThermFEM::SimulationProperties::Instance().setCalculationParameters(
      true, false, false, false, false);
    HygroThermFEM::MultiDomain multiDomain({.performThermal = false, .performMoisture = true});
    buildSandstoneStrip(multiDomain);

    EXPECT_THROW(static_cast<void>(multiDomain.steadyState(HygroThermFEM::SteadyMethod::direct)),
                 std::runtime_error);

    HygroThermFEM::SimulationProperties::Instance().resetCalculationParameters();
}

//! The same refusal when the thermal domain beside it is well posed: prescribed temperatures
//! on both faces and no humidity anywhere. The thermal field has its answer, the moisture
//! field has none, and the march still solves the pair from its start.
TEST(SteadyStateDirect, SealedMoistureBesideThermalIsRefusedAndTheMarchIsNot)
{
    HygroThermFEM::SimulationProperties::Instance().setCalculationParameters(
      true, true, true, true, false);
    HygroThermFEM::MultiDomain multiDomain;
    auto builder = buildSandstoneStrip(multiDomain);
    builder.applyBC_FixedTemperature(TestHelper::BeamBuilder::Edge::Left, -10.0);
    builder.applyBC_FixedTemperature(TestHelper::BeamBuilder::Edge::Right, 20.0);

    EXPECT_TRUE(multiDomain.thermal().hasBoundaryExchange());
    EXPECT_FALSE(multiDomain.moisture().hasBoundaryExchange());
    EXPECT_THROW(static_cast<void>(multiDomain.steadyState(HygroThermFEM::SteadyMethod::direct)),
                 std::runtime_error);

    const auto march = multiDomain.steadyState(HygroThermFEM::SteadyMethod::march);
    EXPECT_TRUE(march.converged);

    HygroThermFEM::SimulationProperties::Instance().resetCalculationParameters();
}
