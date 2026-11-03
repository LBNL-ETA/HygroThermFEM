#include <algorithm>
#include <cmath>
#include <ranges>
#include <vector>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

#include "BeamBuilder.hxx"
#include "TestHelpers.hxx"
#include "TestMaterials.hxx"

// ThermSample_* family: wall assemblies of the kind THERM users build, reproduced as direct
// HygroThermFEM runs. Reproducing them here rather than driving THERM keeps the engine's
// behaviour on realistic geometry under test with no mesher, file format or GUI in the way,
// and makes each case self-contained: the geometry is built in code and the materials come
// from the shared test library (TestMaterials.hxx), so nothing is read at run time.
//
// This case is a two-layer wall, reduced to the 1D beam it effectively is: Stucco 76.3 mm on
// the exterior side, Fiberglass Batts 80 mm on the interior side, 72.8 mm tall, one element
// high. The films are the convection parts of the NFRC 100-2010 exterior and interior
// boundary conditions, held constant: -18 C / RH 0.5 / hc 26 W/m2K on the stucco face and
// 21 C / RH 0.5 / hc 4.65 W/m2K on the fiberglass face. A THERM exterior boundary condition
// adds black-body radiation on top of the film, which is not reproduced here, so absolute
// levels differ a little from an equivalent GUI run; every property checked below is
// internal to this file and unaffected by that difference.
//
// Why this assembly: the cold stucco drives the wall into condensation, so a saturated band
// forms and the humidity bound RH <= 1 becomes active over part of the mesh. That is the
// regime in which the moisture solver is hardest, and in which it has been wrong before.
//
// The criterion is that under constant boundary conditions the steady solution and a
// transient run that has come to rest must be the same field: both discretize the same
// equations and differ only by a time derivative that vanishes at rest. Two forms of that
// statement are used here. The FIXED-POINT form feeds the steady field to the transient
// solver, which must survive an equilibrium-scale timestep unchanged. The REACHABILITY form
// marches the transient from the initial condition until it arrives at the steady field;
// that is the stronger statement, and the one a solver which accepts partial steps fails.
// Both are checked below. Reachability is slow in wall-clock terms only because the approach
// itself is slow in physical terms: the mild isothermal case gets there in a few steps,
// while the saturated case needs about eleven simulated years, the cold stucco storing a
// large amount of moisture that a small vapour flux has to fill.
//
// The same comparison, repeated for each moisture modelling option separately, lives in
// SteadyVsTransient_StuccoFiberglassBeam.
namespace
{
    struct WallFields
    {
        std::vector<double> temperature;
        std::vector<double> humidity;
        std::vector<double> waterContent;
    };

    HygroThermFEM::State wallInitialState()
    {
        // ConstantInitialConditionsSteadyState in the THERM file: 21 C / RH 0.5.
        return {.temperature = 21.0, .humidity = 0.5, .pressure = 101325.0, .liquidPercent = 1.0};
    }

    void buildWall(HygroThermFEM::MultiDomain & multiDomain,
                   const HygroThermFEM::FixedBCHCCoefficients & exterior,
                   const HygroThermFEM::FixedBCHCCoefficients & interior)
    {
        const auto & stucco = multiDomain.materials().createSolidMaterial(TestHelper::Stucco());
        const auto & fiberglass =
          multiDomain.materials().createSolidMaterial(TestHelper::FiberglassBatts());

        TestHelper::BeamBuilder builder(multiDomain);
        builder.xStart(0.0)
          .height(0.0728)
          .numElementsY(1)
          .state(wallInitialState())
          .addSegment({.material = stucco.name(), .numElementsX = 8, .width = 0.0763})
          .addSegment({.material = fiberglass.name(), .numElementsX = 8, .width = 0.08})
          .build();

        builder.applyBC_FixedHc(TestHelper::BeamBuilder::Edge::Left, exterior);
        builder.applyBC_FixedHc(TestHelper::BeamBuilder::Edge::Right, interior);
    }

    //! The THERM file's NFRC films: cold-side saturation band, exercises the bound pinning.
    void buildNfrcWall(HygroThermFEM::MultiDomain & multiDomain)
    {
        buildWall(multiDomain, {-18.0, 26.0, 0.5}, {21.0, 4.65, 0.5});
    }

    //! Mild isothermal case: RH 0.10 / 0.70 at 21 C on both faces. Pure vapor diffusion, no
    //! saturation and no pinning anywhere -- the plain equivalence is provable here before
    //! the NFRC case adds the bounded condensation zone.
    void buildMildWall(HygroThermFEM::MultiDomain & multiDomain)
    {
        buildWall(multiDomain, {21.0, 26.0, 0.10}, {21.0, 4.65, 0.70});
    }

    //! Mirrors ThermDoc::applyEngineSimulationFlags for a steady-state moisture run with the
    //! "water liquid transportation" checkbox off: the three moisture-to-heat coupling terms
    //! are always excluded in steady mode, and liquid transport is excluded here because its
    //! redistribution creeps for simulated years -- vapor-only keeps the fixed point sharp.
    void applyThermSteadyFlags()
    {
        HygroThermFEM::SimulationProperties::Instance().setCalculationParameters(
          true, true, true, true, false);
    }

    double maxAbsDifference(const std::vector<double> & lhs, const std::vector<double> & rhs)
    {
        double result{0.0};
        for(std::size_t index = 0u; index < lhs.size(); ++index)
        {
            result = (std::max)(result, std::abs(lhs[index] - rhs[index]));
        }
        return result;
    }

    //! Every second node: the wall is one element tall, so nodes 2k and 2k+1 share a column.
    //! Fields are compared per column, which is what the recorded tables below hold.
    std::vector<double> perColumn(const std::vector<double> & nodal)
    {
        std::vector<double> result;
        for(std::size_t index = 0u; index < nodal.size(); index += 2u)
        {
            result.push_back(nodal[index]);
        }
        return result;
    }

    struct SteadyRun
    {
        WallFields fields;
        bool converged{false};
        std::size_t passes{0u};
    };

    template<typename BuildWall>
    SteadyRun runSteadyState(const BuildWall & build)
    {
        HygroThermFEM::MultiDomain multiDomain;
        build(multiDomain);
        const auto solution = multiDomain.steadyState();
        return {{solution.temperature, solution.humidity, solution.waterContent},
                solution.converged,
                solution.steadyPasses};
    }

    //! A steady solve that stopped on its pass budget has not solved anything, so every test
    //! that uses a steady field checks this before comparing numbers.
    void expectSteadyConverged(const SteadyRun & run)
    {
        EXPECT_TRUE(run.converged)
          << "steady solve stopped on its pass budget after " << run.passes << " passes";
    }

    void expectMatchesRecorded(const std::vector<double> & actual,
                               const std::vector<double> & recorded,
                               const double tolerance,
                               const char * what)
    {
        ASSERT_EQ(actual.size(), recorded.size()) << what;
        for(std::size_t column = 0u; column < recorded.size(); ++column)
        {
            EXPECT_NEAR(actual[column], recorded[column], tolerance)
              << what << ", column " << column + 1u << " of " << recorded.size();
        }
    }

    // -----------------------------------------------------------------------------------
    // Recorded fields, per x-column from the exterior face (column 1) to the interior face
    // (column 17); column 9 is the stucco/fiberglass interface. Regenerate by running these
    // tests with the environment variable HTF_DUMP_GOLDEN set, which makes the dumpGolden
    // calls below print each array as a C++ initializer, and pasting the result back here.
    // -----------------------------------------------------------------------------------

    //! Mild isothermal wall, 21 C on both faces: a smooth vapour-diffusion profile between
    //! the two film humidities, no saturation anywhere.
    const std::vector<double> mildSteadyHumidity{0.100482187945, 0.170415297324, 0.240348406704, 0.310281516083, 0.380214625462, 0.450147734842, 0.520080844221, 0.590013953601, 0.65994706298, 0.664616667028, 0.669286271075, 0.673955875123, 0.678625479171, 0.683295083218, 0.687964687266, 0.692634291314, 0.697303895361};

    //! NFRC wall. The same model, options and fields as the "all options off" case of
    //! SteadyVsTransient_StuccoFiberglassBeam, which checks them against a ten-year
    //! transient run; here they are the regression guard on the steady solve itself.
    const std::vector<double> nfrcSteadyTemperature{-17.3193405246, -17.1207681324, -16.9221957401, -16.7236233479, -16.5250509557, -16.3264785634, -16.1279061712, -15.9293337789, -15.7307613867, -11.6151459542, -7.49953052172, -3.38391508924, 0.731700343246, 4.84731577573, 8.96293120822, 13.0785466407, 17.1941620732};
    const std::vector<double> nfrcSteadyHumidity{0.473781878417, 0.548040124616, 0.619738081018, 0.688946630249, 0.755734769933, 0.820169663192, 0.882316687755, 0.942239483718, 1, 1, 1, 1, 0.949135976199, 0.86604353093, 0.771978425694, 0.678306528699, 0.590715597351};
}   // namespace

TEST(ThermSample_StuccoFiberglass, MildGradientTransientReachesSteady)
{
    // Tier one: the isothermal 10% / 70% case. With no saturation involved the transient
    // marched from the uniform 50% initial condition must land on the steady solve, which is
    // the reachability statement in its simplest setting, and the steady field itself must
    // match the recorded profile.
    constexpr double dTime{3.6e6};
    constexpr std::size_t maxSteps{50u};
    constexpr double agreementTolerance{1e-3};
    constexpr double recordedTolerance{1e-6};

    applyThermSteadyFlags();
    const auto steady = runSteadyState(buildMildWall);

    HygroThermFEM::MultiDomain multiDomain;
    buildMildWall(multiDomain);
    auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
    auto humidities = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);

    std::size_t stepsTaken{maxSteps};
    for(std::size_t step = 0u; step < maxSteps; ++step)
    {
        const auto solution = multiDomain.transient(temperatures, humidities, dTime, step);
        temperatures = solution.temperature;
        humidities = solution.humidity;
        if(maxAbsDifference(humidities, steady.fields.humidity) < agreementTolerance)
        {
            stepsTaken = step + 1u;
            break;
        }
    }
    HygroThermFEM::SimulationProperties::Instance().resetCalculationParameters();

    expectSteadyConverged(steady);
    TestHelper::dumpGolden("mildSteadyHumidity", perColumn(steady.fields.humidity));
    expectMatchesRecorded(perColumn(steady.fields.humidity),
                          mildSteadyHumidity,
                          recordedTolerance,
                          "mild wall steady humidity");

    EXPECT_LT(stepsTaken, maxSteps)
      << "transient did not reach the steady field within " << maxSteps << " steps; it ended "
      << maxAbsDifference(humidities, steady.fields.humidity) << " RH away";
    EXPECT_LT(maxAbsDifference(temperatures, steady.fields.temperature), agreementTolerance);
    for(const auto value : humidities)
    {
        EXPECT_GE(value, 0.0);
        EXPECT_LE(value, 1.0);
    }
}

TEST(ThermSample_StuccoFiberglass, SteadyStateBoundedProfile)
{
    // The bounded steady solve on the real wall, against the recorded profile. The two
    // properties that give the numbers their meaning are asserted separately, so a failure
    // says which one broke: humidity must respect [0, 1] (the unconstrained Glaser solve
    // supersaturates half the wall) and a saturated condensation band must form at the cold
    // end of the vapour-open fiberglass.
    constexpr double recordedTolerance{1e-6};

    applyThermSteadyFlags();
    const auto steady = runSteadyState(buildNfrcWall);
    HygroThermFEM::SimulationProperties::Instance().resetCalculationParameters();

    expectSteadyConverged(steady);
    TestHelper::dumpGolden("nfrcSteadyTemperature", perColumn(steady.fields.temperature));
    TestHelper::dumpGolden("nfrcSteadyHumidity", perColumn(steady.fields.humidity));
    expectMatchesRecorded(perColumn(steady.fields.temperature),
                          nfrcSteadyTemperature,
                          recordedTolerance,
                          "NFRC wall steady temperature");
    expectMatchesRecorded(perColumn(steady.fields.humidity),
                          nfrcSteadyHumidity,
                          recordedTolerance,
                          "NFRC wall steady humidity");

    for(const auto value : steady.fields.humidity)
    {
        EXPECT_GE(value, 0.0);
        EXPECT_LE(value, 1.0 + 1e-9);
    }
    EXPECT_TRUE(std::ranges::any_of(steady.fields.humidity,
                                    [](const double value) { return value > 0.999; }))
      << "expected a saturated condensation band";
}

TEST(ThermSample_StuccoFiberglass, SteadyStateIsTransientFixedPoint)
{
    // The documented equivalence: seed the transient solver with the steady fields and take
    // one equilibrium-scale step -- the state must not move. Any disagreement between the
    // two solvers' discretizations would push the solution away from the seed.
    constexpr double dTime{3.6e6};
    // The step leaves the field where it found it: any residual is the two solvers' own
    // iteration tolerances, not a difference in what they discretize.
    constexpr double agreementTolerance{5e-3};

    applyThermSteadyFlags();
    const auto steady = runSteadyState(buildNfrcWall);

    HygroThermFEM::MultiDomain multiDomain;
    buildNfrcWall(multiDomain);
    multiDomain.nodes().updateNodeTemperatures(steady.fields.temperature, true);
    multiDomain.nodes().updateNodeHumidities(steady.fields.humidity, true);
    const auto solution =
      multiDomain.transient(steady.fields.temperature, steady.fields.humidity, dTime, 0u);
    HygroThermFEM::SimulationProperties::Instance().resetCalculationParameters();

    expectSteadyConverged(steady);
    const double temperatureDrift{
      maxAbsDifference(solution.temperature, steady.fields.temperature)};
    const double humidityDrift{maxAbsDifference(solution.humidity, steady.fields.humidity)};

    EXPECT_LT(temperatureDrift, agreementTolerance)
      << "steady field moved over one " << dTime / 3600.0 << " h step";
    EXPECT_LT(humidityDrift, agreementTolerance)
      << "steady field moved over one " << dTime / 3600.0 << " h step";
}

TEST(ThermSample_StuccoFiberglass, SaturatedTransientReachesSteady)
{
    // Tier two: reachability on the wall's own saturated case, where a condensation band
    // forms and the humidity bound is active across part of the mesh. The transient is
    // marched from the uniform initial condition with equilibrium-scale steps and must
    // arrive at the steady field.
    //
    // This is the test that a solver accepting partial steps fails, and it did fail until
    // the moisture Newton took a true Newton direction with the bound applied to the
    // correction: the approach used to stall at a per-step creep that neither a larger
    // timestep nor a larger iteration budget could remove, ending 0.23 RH away. It now
    // contracts geometrically, roughly halving every ten steps.
    constexpr double dTime{3.6e6};   // 1000 h; the 100 steps below are ~11 simulated years
    constexpr std::size_t maxSteps{100u};
    constexpr double humidityTolerance{1e-3};
    constexpr double temperatureTolerance{1e-2};

    applyThermSteadyFlags();
    const auto steady = runSteadyState(buildNfrcWall);

    HygroThermFEM::MultiDomain multiDomain;
    buildNfrcWall(multiDomain);
    auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
    auto humidities = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);
    for(std::size_t step = 0u; step < maxSteps; ++step)
    {
        const auto solution = multiDomain.transient(temperatures, humidities, dTime, step);
        temperatures = solution.temperature;
        humidities = solution.humidity;
    }
    HygroThermFEM::SimulationProperties::Instance().resetCalculationParameters();

    expectSteadyConverged(steady);
    EXPECT_LT(maxAbsDifference(humidities, steady.fields.humidity), humidityTolerance)
      << "after " << maxSteps << " x " << dTime / 3600.0 << " h steps";
    EXPECT_LT(maxAbsDifference(temperatures, steady.fields.temperature), temperatureTolerance)
      << "after " << maxSteps << " x " << dTime / 3600.0 << " h steps";

    // The endpoint is the recorded steady profile itself, so this test also pins the field
    // the transient arrives at, not only its distance from a value computed in the same run.
    expectMatchesRecorded(perColumn(humidities),
                          nfrcSteadyHumidity,
                          humidityTolerance,
                          "NFRC wall transient endpoint humidity");
}
