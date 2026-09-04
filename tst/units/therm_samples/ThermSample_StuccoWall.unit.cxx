#include <algorithm>
#include <array>
#include <limits>
#include <vector>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

#include "BeamBuilder.hxx"
#include "TestHelpers.hxx"
#include "TestMaterials.hxx"

// ThermSample_* family: real-world THERM sample files reproduced as direct HygroThermFEM
// runs, to be exercised (gtest_filter=ThermSample_*) BEFORE any THERM build consumes a new
// HygroThermFEM. Source of truth for the geometry/BC values: the .thmz sample noted below.
//
// The THERM sample "Stucco Wall - Moisture" (Samples\Stucco Wall - Moisture.thmz) as a
// 1D beam: Stucco 25.4 mm | Laminated panel 19.05 mm | Fiberglass Batts 76.2 mm |
// Gypsum Board Interior 25.4 mm, exterior film 10 C / RH 0.10 / hc 7 on the stucco side,
// interior film 40 C / RH 0.30 / hc 10 on the gypsum side ("use first timestep" BCs),
// dt = 360 s. This is the real-world case reported to blow up when started from a
// near-saturated initial humidity (0.99); the mild sample initial condition (0.10) is
// kept as the companion case.
namespace
{
    struct WallRunResult
    {
        bool completed = false;
        std::string error;
        std::vector<std::vector<double>> temperatures;
        std::vector<std::vector<double>> waterContents;
        std::vector<std::vector<double>> humidities;
    };

    WallRunResult runStuccoWall(const double initialHumidity,
                                const unsigned numberOfSteps,
                                const double dTime = 360.0,
                                const unsigned meshScale = 1,
                                const std::array<unsigned, 4> & elementsPerLayer = {4, 3, 6, 4},
                                const unsigned numElementsY = 1)
    {
        WallRunResult result;
        HygroThermFEM::MultiDomain multiDomain;

        const HygroThermFEM::State initialState({.temperature = 21.0,
                                                 .humidity = initialHumidity,
                                                 .pressure = 101325.0,
                                                 .liquidPercent = 1.0});

        const auto & stucco = multiDomain.materials().createSolidMaterial(TestHelper::Stucco());
        const auto & panel =
          multiDomain.materials().createSolidMaterial(TestHelper::LaminatedPanel());
        const auto & fiberglass =
          multiDomain.materials().createSolidMaterial(TestHelper::FiberglassBatts());
        const auto & gypsum =
          multiDomain.materials().createSolidMaterial(TestHelper::GypsumBoardInterior());

        TestHelper::BeamBuilder builder(multiDomain);
        builder.xStart(0.0)
          .height(0.05)
          .numElementsY(numElementsY)
          .state(initialState)
          .addSegment({.material = stucco.name(),
                       .numElementsX = elementsPerLayer[0] * meshScale,
                       .width = 0.0254})
          .addSegment({.material = panel.name(),
                       .numElementsX = elementsPerLayer[1] * meshScale,
                       .width = 0.01905})
          .addSegment({.material = fiberglass.name(),
                       .numElementsX = elementsPerLayer[2] * meshScale,
                       .width = 0.0762})
          .addSegment({.material = gypsum.name(),
                       .numElementsX = elementsPerLayer[3] * meshScale,
                       .width = 0.0254})
          .build();

        const HygroThermFEM::FixedBCHCCoefficients exterior{10.0, 7.0, 0.1};
        const HygroThermFEM::FixedBCHCCoefficients interior{40.0, 10.0, 0.3};
        builder.applyBC_FixedHc(TestHelper::BeamBuilder::Edge::Left, exterior);
        builder.applyBC_FixedHc(TestHelper::BeamBuilder::Edge::Right, interior);

        try
        {
            auto temperatures =
              multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
            auto humidities = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);

            result.temperatures.push_back(temperatures);
            result.humidities.push_back(humidities);
            result.waterContents.push_back(
              multiDomain.nodes().properties(HygroThermFEM::Variable::water));

            for(unsigned step = 0; step < numberOfSteps; ++step)
            {
                const auto solution =
                  multiDomain.transient(temperatures, humidities, dTime, step);
                temperatures = solution.temperature;
                humidities = solution.humidity;
                result.temperatures.push_back(temperatures);
                result.humidities.push_back(humidities);
                result.waterContents.push_back(solution.waterContent);
            }
            result.completed = true;
        }
        catch(const std::exception & exc)
        {
            result.error = exc.what();
        }

        return result;
    }

    struct OscillationSummary
    {
        std::size_t maxExtrema{0};
        double maxAmplitude{0.0};
        std::size_t worstStep{0};
    };

    //! Detects spatial sawtooth oscillations in a per-step nodal field. Extracts the
    //! bottom node row (column-major layout, stride = numElementsY + 1) as the x-ordered
    //! profile and flags interior local extrema; the amplitude of an extremum is the
    //! smaller jump to its neighbours, so a single-node kink registers at its full size.
    OscillationSummary analyzeOscillations(const std::vector<std::vector<double>> & fields,
                                           const std::size_t stride = 2)
    {
        OscillationSummary summary;
        for(std::size_t step = 0; step < fields.size(); ++step)
        {
            std::vector<double> profile;
            for(std::size_t idx = 0; idx < fields[step].size(); idx += stride)
            {
                profile.push_back(fields[step][idx]);
            }
            std::size_t extrema{0};
            double amplitude{0.0};
            for(std::size_t idx = 1; idx + 1 < profile.size(); ++idx)
            {
                const double diffLeft = profile[idx] - profile[idx - 1];
                const double diffRight = profile[idx + 1] - profile[idx];
                if(diffLeft * diffRight < 0.0)
                {
                    ++extrema;
                    amplitude = (std::max)(
                      amplitude, (std::min)(std::abs(diffLeft), std::abs(diffRight)));
                }
            }
            if(amplitude > summary.maxAmplitude)
            {
                summary.maxAmplitude = amplitude;
                summary.worstStep = step;
            }
            summary.maxExtrema = (std::max)(summary.maxExtrema, extrema);
        }
        return summary;
    }

    //! Where a run reaches its extremes, so a breach can be reported by step and node.
    struct SeriesExtremes
    {
        double minimum{(std::numeric_limits<double>::max)()};
        double maximum{std::numeric_limits<double>::lowest()};
        std::size_t minimumStep{0u};
        std::size_t minimumNode{0u};
        std::size_t maximumStep{0u};
        std::size_t maximumNode{0u};
    };

    SeriesExtremes extremesOf(const std::vector<std::vector<double>> & series)
    {
        SeriesExtremes extremes;
        for(std::size_t step = 0u; step < series.size(); ++step)
        {
            for(std::size_t node = 0u; node < series[step].size(); ++node)
            {
                const double value{series[step][node]};
                if(value < extremes.minimum)
                {
                    extremes.minimum = value;
                    extremes.minimumStep = step;
                    extremes.minimumNode = node;
                }
                if(value > extremes.maximum)
                {
                    extremes.maximum = value;
                    extremes.maximumStep = step;
                    extremes.maximumNode = node;
                }
            }
        }
        return extremes;
    }

    //! Two assertions over the whole run rather than one per value: a breach then reports
    //! where it happened and how far out it went, which is what a blow-up needs, and a run
    //! that stays inside says nothing at all.
    void expectWithinEnvelope(const std::vector<std::vector<double>> & series,
                              const double low,
                              const double high,
                              const char * what)
    {
        const auto extremes = extremesOf(series);
        EXPECT_GE(extremes.minimum, low)
          << what << " fell below the envelope at step " << extremes.minimumStep << ", node "
          << extremes.minimumNode;
        EXPECT_LE(extremes.maximum, high)
          << what << " rose above the envelope at step " << extremes.maximumStep << ", node "
          << extremes.maximumNode;
    }
}   // namespace

TEST(ThermSample_StuccoWall, SampleInitialConditions)
{
    // The stock sample: initial 21 C / RH 0.10. Must run to completion with bounded,
    // physical fields.
    const auto result = runStuccoWall(0.1, 100);
    ASSERT_TRUE(result.completed) << result.error;

    // The two film temperatures plus a transient margin.
    expectWithinEnvelope(result.temperatures, 10.0 - 2.0, 40.0 + 2.0, "temperature");
    expectWithinEnvelope(result.humidities, 0.0, 1.0, "humidity");
}

TEST(ThermSample_StuccoWall, NearSaturatedNoLatentHeat)
{
    // Discriminator for the latent-heat sign defect: the same near-saturated start with
    // the heat-of-evaporation term excluded. If this stays bounded while the full-physics
    // variant self-heats above both boundary temperatures, the interior latent term is
    // the driver (heatOfEvaporation() dropped the minus sign of Theoretical Model eq. 59).
    HygroThermFEM::SimulationProperties::Instance().setCalculationParameters(
      false, true, false, false, false);
    const auto result = runStuccoWall(0.99, 100);
    HygroThermFEM::SimulationProperties::Instance().resetCalculationParameters();

    ASSERT_TRUE(result.completed) << result.error;
    expectWithinEnvelope(
      result.temperatures, 10.0 - 2.0, 40.0 + 2.0, "temperature without the latent term");
}

TEST(ThermSample_StuccoWall, NearSaturatedFineMesh)
{
    // The "Linear solver failed to factorize the system matrix" reported from THERM: this
    // case, on an 8x finer mesh (~136 elements), started from humidities right up against
    // saturation. It was a live defect and this was a printing probe for it.
    //
    // It stopped reproducing once the moisture Newton took a true Newton direction with the
    // humidity bound applied to the correction system (2026-09-03), so the case is kept as a
    // guard: all three starts must now run to completion with a physical field.
    for(const double phi0 : {0.99, 0.999, 0.99999})
    {
        const auto result = runStuccoWall(phi0, 20, 360.0, 8);
        ASSERT_TRUE(result.completed) << "initial humidity " << phi0 << ": " << result.error;
        expectWithinEnvelope(result.temperatures, 10.0, 40.0, "temperature");
        expectWithinEnvelope(result.humidities, 0.0, 1.0, "humidity");
    }
}

TEST(ThermSample_StuccoWall, SaturatedThermMeshNoSawtooth)
{
    // Regression guard for the explicit-latent-term sawtooth (fixed by the chain-rule
    // split: T-part implicit in the conductance). The exact THERM quadtree element counts
    // (3/2/7/3) at full saturation were the reported failure: T sawtooth 12.9 K @ step 77
    // before the fix. The full mesh/equation sweep lives in DISABLED_SaturatedMeshBisect.
    const auto result = runStuccoWall(1.0, 100, 360.0, 1, {3, 2, 7, 3});
    ASSERT_TRUE(result.completed) << result.error;

    const auto tempOsc = analyzeOscillations(result.temperatures);
    EXPECT_LT(tempOsc.maxAmplitude, 0.1)
      << "spatial temperature sawtooth at step " << tempOsc.worstStep;
    for(const auto & temps : result.temperatures)
    {
        for(const auto val : temps)
        {
            EXPECT_GT(val, -5.0);
            EXPECT_LT(val, 45.0);
        }
    }
}

TEST(ThermSample_StuccoWall, SaturatedMeshRefinementStaysStable)
{
    // The 2026-07-13 THERM-Viz report: at exact saturation on a fine mesh the fields blew up
    // around step 54 while a coarse mesh stayed bounded. It stopped reproducing once the
    // moisture Newton took a true Newton direction with the humidity bound applied to the
    // correction (2026-09-03), so the mesh sweep is kept as a guard: refining the mesh must
    // not destabilise the run.
    //
    // Its sibling SaturatedThermMeshNoSawtooth covers the mesh THERM actually generated for
    // that report, 3/2/7/3 elements per layer; this one covers uniform refinement, where the
    // blow-up was reported. The per-term sweep that hunted the driver stays disabled below.
    constexpr double phiStart{1.0};
    constexpr unsigned steps{100};

    for(const unsigned meshScale : {1u, 8u})
    {
        const auto result = runStuccoWall(phiStart, steps, 360.0, meshScale);
        ASSERT_TRUE(result.completed) << "mesh x" << meshScale << ": " << result.error;

        // Measured across these scales: 4.1 to 24.8 C. The envelope is the one its sibling
        // uses, wide enough to allow the evaporative cooling this fully saturated start
        // produces and narrow enough to catch the reported blow-up.
        expectWithinEnvelope(result.temperatures, -5.0, 45.0, "temperature");
        expectWithinEnvelope(result.humidities, 0.0, 1.0, "humidity");

        const auto tempOsc = analyzeOscillations(result.temperatures);
        EXPECT_LT(tempOsc.maxAmplitude, 0.1)
          << "mesh x" << meshScale << ": spatial temperature sawtooth at step "
          << tempOsc.worstStep;
    }
}

TEST(ThermSample_StuccoWall, DISABLED_SaturatedMeshBisect)
{
    // Diagnostic kept for the open finding it carries: with heat of evaporation excluded at
    // mesh x8 the wall reaches 47.3 C, above the 40 C interior film, first leaving the
    // envelope at step 50. The vapour and capillary conduction terms still advect enthalpy
    // when the latent term is off, so energy is carried without its latent counterpart --
    // the same shape of problem as capillary conduction running without a liquid flux. It is
    // mesh-dependent: the coarse-mesh version of that case passes. Heat of evaporation is a
    // checkbox a THERM user can untick, so this is reachable.
    //
    // Excludes one term at a time at the fine scale to find the driver. Run explicitly with
    // --gtest_also_run_disabled_tests.
    constexpr double phiStart{1.0};
    constexpr unsigned steps{100};

    const auto report = [](const std::string & label, const WallRunResult & result) {
        double minTemp{1e9};
        double maxTemp{-1e9};
        size_t firstBad{0};
        bool haveBad{false};
        for(size_t step = 0; step < result.temperatures.size(); ++step)
        {
            for(const auto val : result.temperatures[step])
            {
                minTemp = (std::min)(minTemp, val);
                maxTemp = (std::max)(maxTemp, val);
                if(!haveBad && (val < -5.0 || val > 45.0))
                {
                    firstBad = step;
                    haveBad = true;
                }
            }
        }
        const auto tempOsc = analyzeOscillations(result.temperatures);
        const auto phiOsc = analyzeOscillations(result.humidities);
        std::cout << "[SatBisect] " << label
                  << (result.completed ? "" : " ABORTED: " + result.error) << " T=["
                  << minTemp << ", " << maxTemp << "]"
                  << (haveBad ? " first out-of-envelope step " + std::to_string(firstBad)
                              : "")
                  << " | T osc: n=" << tempOsc.maxExtrema << " amp=" << tempOsc.maxAmplitude
                  << " @step " << tempOsc.worstStep
                  << " | phi osc: n=" << phiOsc.maxExtrema << " amp=" << phiOsc.maxAmplitude
                  << " @step " << phiOsc.worstStep << std::endl;
    };

    struct Exclusion
    {
        std::string label;
        bool liquid{false};
        bool evaporation{false};
        bool capillary{false};
        bool vapor{false};
    };
    const std::vector<Exclusion> exclusions{
      {.label = "no liquid transport", .liquid = true},
      {.label = "no heat of evaporation", .evaporation = true},
      {.label = "no capillary conduction", .capillary = true},
      {.label = "no vapor diffusion conduction", .vapor = true}};
    for(const auto & exclusion : exclusions)
    {
        HygroThermFEM::SimulationProperties::Instance().setCalculationParameters(
          exclusion.liquid,
          exclusion.evaporation,
          exclusion.capillary,
          exclusion.vapor,
          false);
        report(exclusion.label + ", mesh x8", runStuccoWall(phiStart, steps, 360.0, 8));
        HygroThermFEM::SimulationProperties::Instance().resetCalculationParameters();
    }
}

TEST(ThermSample_StuccoWall, NearSaturatedRefinedDt)
{
    // The same near-saturated case at a quarter of the timestep, dt = 90 s. This began as a
    // probe asking whether the sub-freezing surface dip and the latent overshoot the wall
    // once showed were artefacts of the coarse step. They were: both are gone, and refining
    // the step no longer changes the answer materially, 10.1 to 35.3 C here against 11.6 to
    // 35.3 at dt = 360 s.
    //
    // Kept as a guard on that. The bound is a degree below the exterior film rather than the
    // film itself, because the measured floor sits only 0.06 C above it and evaporative
    // cooling legitimately pushes that way; it is still far tighter than the 2.2 C dip and
    // 41.3 C overshoot this case used to produce.
    const auto result = runStuccoWall(0.99, 400, 90.0);
    ASSERT_TRUE(result.completed) << result.error;

    expectWithinEnvelope(result.temperatures, 9.0, 40.0, "temperature");
    expectWithinEnvelope(result.humidities, 0.0, 1.0, "humidity");
}

TEST(ThermSample_StuccoWall, NearSaturatedInitialConditions)
{
    // The formerly reported blow-up: the wall started from RH 0.99. Before the latent-term
    // corrections (sign of h_lg per Theoretical Model eq. 59; flux-consistent vapour
    // potential c_sat*phi) this self-heated to the 1000-clamp; now temperatures stay
    // within the two film temperatures: the field now runs between 11.6 and 35.3 C at this
    // dt = 360 s, and between 10.1 and 35.3 C at the finer dt = 90 s of
    // NearSaturatedRefinedDt.
    //
    // The bound below is therefore the driving temperatures themselves rather than the wide
    // envelope this test used to carry. Latent effects can in principle break it in both
    // directions, evaporative cooling pulling the surface toward the exterior air's wet-bulb
    // temperature and condensation warming the interior side above its film, and they used to
    // do exactly that. If this starts failing, that is the signal, and the bound should be
    // reopened with a measurement rather than widened on principle. NOTE: transients crossing
    // 0 C would make the (currently disabled) freezing model (D5) relevant here.
    const auto result = runStuccoWall(0.99, 100);
    ASSERT_TRUE(result.completed) << result.error;

    expectWithinEnvelope(result.temperatures, 10.0, 40.0, "temperature");
    expectWithinEnvelope(result.humidities, 0.0, 1.0, "humidity");
}
