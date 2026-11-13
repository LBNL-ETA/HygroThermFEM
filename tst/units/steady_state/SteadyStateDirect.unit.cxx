#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

#include "BeamBuilder.hxx"
#include "LinearSolver.hxx"
#include "TestMaterials.hxx"

// The steady solve (MultiDomain::steadyState) on a small case with a recorded field, a sealed
// moisture part, a partly sealed one, and the case it must refuse.
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

    //! Water of the sandstone strip per unit height and depth [kg/m]: the nodal water
    //! content summed over the nodal lengths, half an element at each end column. The two
    //! nodes of a column carry the same value, so the bottom one stands for it.
    double stripWater(const std::vector<double> & waterContent)
    {
        constexpr std::size_t nodesPerColumn{2u};
        constexpr double elementWidth{0.1 / 4.0};
        const std::size_t columns{waterContent.size() / nodesPerColumn};
        double water{0.0};
        for(std::size_t column = 0u; column < columns; ++column)
        {
            const bool atEnd{column == 0u || column + 1u == columns};
            water += (atEnd ? 0.5 : 1.0) * elementWidth * waterContent[column * nodesPerColumn];
        }
        return water;
    }
}   // namespace

//! A linear thermal problem: one pass, one linear solve, the recorded field to solver precision.
TEST(SteadyStateDirect, LinearConductionMatchesRecordedField)
{
    HygroThermFEM::MultiDomain multiDomain({.performThermal = true, .performMoisture = false});
    buildConductionStrip(multiDomain);

    const auto solution = multiDomain.steadyState();

    EXPECT_TRUE(solution.converged);
    const std::vector<double> recorded{
      -17.87392241, -17.87392241, 7.341594828, 7.341594828, 19.94935345, 19.94935345};
    ASSERT_EQ(solution.temperature.size(), recorded.size());
    for(std::size_t node = 0u; node < recorded.size(); ++node)
    {
        EXPECT_NEAR(solution.temperature[node], recorded[node], 1e-6) << "node " << node + 1u;
    }
}

//! A coupled solve whose moisture is sealed: prescribed temperatures on both faces, no
//! humidity anywhere. The thermal field has its answer from the faces; the moisture field
//! has its profile from the flux equations and its level from the water of the start.
TEST(SteadyStateDirect, SealedMoistureBesideThermalKeepsItsWater)
{
    HygroThermFEM::SimulationProperties::Instance().setCalculationParameters(
      true, true, true, true, false);
    HygroThermFEM::MultiDomain multiDomain;
    auto builder = buildSandstoneStrip(multiDomain);
    builder.applyBC_FixedTemperature(TestHelper::BeamBuilder::Edge::Left, 10.0);
    builder.applyBC_FixedTemperature(TestHelper::BeamBuilder::Edge::Right, 30.0);
    const double startWater{stripWater(multiDomain.nodes().properties(HygroThermFEM::Variable::water))};

    const auto solution = multiDomain.steadyState();

    HygroThermFEM::SimulationProperties::Instance().resetCalculationParameters();

    EXPECT_TRUE(solution.converged);
    EXPECT_NEAR(solution.temperature.front(), 10.0, 1e-6);
    EXPECT_NEAR(solution.temperature.back(), 30.0, 1e-6);
    // The cold face is the humid one, and the water is the water the strip started with.
    EXPECT_GT(solution.humidity.front(), solution.humidity.back() + 0.05);
    EXPECT_NEAR(stripWater(solution.waterContent) / startWater, 1.0, 1e-10);
}

//! A strip with a film on one face only is not sealed: the film fixes the level, and the
//! strip goes to the humidity of the air whatever water it started with.
TEST(SteadyStateDirect, PartlySealedMoistureFollowsItsBoundary)
{
    HygroThermFEM::SimulationProperties::Instance().setCalculationParameters(
      true, false, false, false, false);
    HygroThermFEM::MultiDomain multiDomain({.performThermal = false, .performMoisture = true});
    auto builder = buildSandstoneStrip(multiDomain);
    builder.applyBC_FixedHc(TestHelper::BeamBuilder::Edge::Left,
                            HygroThermFEM::FixedBCHCCoefficients{20.0, 5.0, 0.8});
    EXPECT_TRUE(multiDomain.moisture().sealedRegions().empty());

    const auto solution = multiDomain.steadyState();

    HygroThermFEM::SimulationProperties::Instance().resetCalculationParameters();

    EXPECT_TRUE(solution.converged);
    for(const double humidity : solution.humidity)
    {
        EXPECT_NEAR(humidity, 0.8, 1e-6);
    }
}

//! A sealed part whose material stores no water (a sorption curve that is zero at every
//! humidity) has nothing to keep: the solve says so.
TEST(SteadyStateDirect, SealedMoistureWithoutStorageIsRefused)
{
    HygroThermFEM::SimulationProperties::Instance().setCalculationParameters(
      true, false, false, false, false);
    HygroThermFEM::MultiDomain multiDomain({.performThermal = false, .performMoisture = true});
    auto material = TestHelper::LinearSorption();
    material.name = "No storage";
    material.sorptionCurve = {{0.0, 0.0}, {1.0, 0.0}};
    TestHelper::BeamBuilder builder(multiDomain);
    builder.xStart(0.0)
      .height(0.01)
      .numElementsY(1)
      .state({.temperature = 20.0, .humidity = 0.5, .pressure = 101325.0, .liquidPercent = 1.0})
      .addSegment({.material = multiDomain.materials().createSolidMaterial(material).name(),
                   .numElementsX = 4,
                   .width = 0.1})
      .build();

    EXPECT_THROW(static_cast<void>(multiDomain.steadyState()), std::runtime_error);

    HygroThermFEM::SimulationProperties::Instance().resetCalculationParameters();
}

//! A transient step that starts on the steady state has nothing to solve. With the latent
//! heat of fusion modelled, the thermal Newton stops on the reduction of its residual; on
//! this start the residual is already the rounding of the prescribed-temperature rows and
//! cannot be reduced, so the iteration must recognise it instead of running its budget
//! (2026-11-13: about 400 linear solves per step before, on a strip like this one).
TEST(SteadyStateDirect, ATransientStepFromTheSteadyStateTakesFewSolves)
{
    HygroThermFEM::MultiDomain multiDomain;
    multiDomain.setPhysicsOptions({.excludeLatentHeatOfFusion = false});
    auto builder = buildSandstoneStrip(multiDomain);
    builder.applyBC_FixedTemperature(TestHelper::BeamBuilder::Edge::Left, 10.0);
    builder.applyBC_FixedTemperature(TestHelper::BeamBuilder::Edge::Right, 30.0);
    const auto steady = multiDomain.steadyState();
    ASSERT_TRUE(steady.converged);

    constexpr double tenHours{36000.0};
    const auto solvesBefore = HygroThermFEM::CLinearSolver::solveCount();
    const auto step = multiDomain.transient(steady.temperature, steady.humidity, tenHours);
    const auto solves = HygroThermFEM::CLinearSolver::solveCount() - solvesBefore;

    EXPECT_LT(solves, 30u);
    for(std::size_t node = 0u; node < steady.temperature.size(); ++node)
    {
        EXPECT_NEAR(step.temperature[node], steady.temperature[node], 1e-9) << "node " << node + 1u;
        EXPECT_NEAR(step.humidity[node], steady.humidity[node], 1e-9) << "node " << node + 1u;
    }
}
