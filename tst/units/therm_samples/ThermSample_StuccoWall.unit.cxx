#include <algorithm>
#include <array>
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

    double maxAbsValue(const std::vector<std::vector<double>> & values)
    {
        double result = 0.0;
        for(const auto & row : values)
        {
            for(const auto val : row)
            {
                result = (std::max)(result, std::abs(val));
            }
        }
        return result;
    }
}   // namespace

TEST(ThermSample_StuccoWall, SampleInitialConditions)
{
    // The stock sample: initial 21 C / RH 0.10. Must run to completion with bounded,
    // physical fields.
    const auto result = runStuccoWall(0.1, 100);
    ASSERT_TRUE(result.completed) << result.error;

    for(const auto & temps : result.temperatures)
    {
        for(const auto val : temps)
        {
            EXPECT_GT(val, 10.0 - 2.0);
            EXPECT_LT(val, 40.0 + 2.0);
        }
    }
    for(const auto & phis : result.humidities)
    {
        for(const auto val : phis)
        {
            EXPECT_GE(val, 0.0);
            EXPECT_LE(val, 1.0);
        }
    }
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
    std::cout << "[StuccoWall 0.99, no latent] max|T| over run = "
              << maxAbsValue(result.temperatures) << "\n";
    for(const auto & temps : result.temperatures)
    {
        for(const auto val : temps)
        {
            EXPECT_GT(val, 10.0 - 2.0);
            EXPECT_LT(val, 40.0 + 2.0);
        }
    }
}

TEST(ThermSample_StuccoWall, DISABLED_NearSaturatedFineMesh)
{
    // THERM-mesh-density probe: same case, 8x finer mesh (~136 elements). Chases the
    // "Linear solver failed to factorize the system matrix" reported from THERM.
    for(const double phi0 : {0.99, 0.999, 0.99999})
    {
        const auto result = runStuccoWall(phi0, 20, 360.0, 8);
        std::cout << "[FineMesh] phi0=" << phi0
                  << (result.completed ? " completed" : (" FAILED: " + result.error))
                  << (result.completed
                        ? ", max|T|=" + std::to_string(maxAbsValue(result.temperatures))
                        : "")
                  << std::endl;
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

TEST(ThermSample_StuccoWall, DISABLED_SaturatedMeshBisect)
{
    // Diagnostic for the 2026-07-13 THERM-Viz report: phi0 = 1.0 (exact saturation),
    // fine mesh, dt = 360 s -- fields blow up around step 54 while a coarse mesh stays
    // bounded. Stage 1 sweeps mesh scales under full physics to locate the failing
    // scale; stage 2 excludes one term at a time at the fine scale to find the driver.
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

    // The mesh THERM generated for the reported run: 3/2/7/3 elements per layer -- a
    // coarse panel layer against a fine fiberglass layer. Reproduced exactly, alongside
    // the uniform-scale sweep.
    const std::array<unsigned, 4> thermMesh{3, 2, 7, 3};
    report("full physics, THERM mesh 3/2/7/3",
           runStuccoWall(phiStart, steps, 360.0, 1, thermMesh));

    for(const unsigned meshScale : {1u, 2u, 4u, 8u})
    {
        report("full physics, mesh x" + std::to_string(meshScale),
               runStuccoWall(phiStart, steps, 360.0, meshScale));
    }

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

TEST(ThermSample_StuccoWall, DISABLED_NearSaturatedRefinedDt)
{
    // dt-refinement probe (diagnostic, run explicitly with --gtest_also_run_disabled_tests):
    // same case at dt = 90 s. If the sub-freezing surface dip and the +latent overshoot
    // shrink toward the wet-bulb/condensation bounds, they are big-dt transients.
    const auto result = runStuccoWall(0.99, 400, 90.0);
    ASSERT_TRUE(result.completed) << result.error;
    double minTemp = 1e9;
    double maxTemp = -1e9;
    for(const auto & temps : result.temperatures)
    {
        for(const auto val : temps)
        {
            minTemp = (std::min)(minTemp, val);
            maxTemp = (std::max)(maxTemp, val);
        }
    }
    std::cout << "[StuccoWall 0.99, dt=90] T range over run = [" << minTemp << ", "
              << maxTemp << "]" << std::endl;
}

TEST(ThermSample_StuccoWall, NearSaturatedInitialConditions)
{
    // The formerly reported blow-up: the wall started from RH 0.99. Before the latent-term
    // corrections (sign of h_lg per Theoretical Model eq. 59; flux-consistent vapour
    // potential c_sat*phi) this self-heated to the 1000-clamp; now temperatures stay
    // within the latent-effect envelope: evaporative cooling can reach the exterior air's
    // wet-bulb temperature (~2 C for 10 C / RH 0.10, briefly overshot at this coarse
    // dt = 360 s) and condensation warms a few degrees above the warm boundary. The
    // DISABLED_NearSaturatedRefinedDt probe documents convergence to [2.2, 41.3] C at
    // dt = 90 s. NOTE: transients crossing 0 C make the (currently disabled) freezing
    // model (D5) relevant for this wall in colder climates.
    const auto result = runStuccoWall(0.99, 100);
    ASSERT_TRUE(result.completed) << result.error;

    std::cout << "[StuccoWall 0.99] max|T| over run = " << maxAbsValue(result.temperatures)
              << ", final T range = ["
              << *std::min_element(result.temperatures.back().begin(),
                                   result.temperatures.back().end())
              << ", "
              << *std::max_element(result.temperatures.back().begin(),
                                   result.temperatures.back().end())
              << "]\n";

    for(const auto & temps : result.temperatures)
    {
        for(const auto val : temps)
        {
            EXPECT_GT(val, -5.0);   // wet-bulb floor plus coarse-dt transient margin
            EXPECT_LT(val, 45.0);   // condensation-warming ceiling plus margin
        }
    }
    for(const auto & phis : result.humidities)
    {
        for(const auto val : phis)
        {
            EXPECT_GE(val, 0.0);
            EXPECT_LE(val, 1.0);
        }
    }
}
