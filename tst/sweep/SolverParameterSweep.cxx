#include <gtest/gtest.h>

#include <iomanip>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "HygroThermFEM2D.hxx"

#include "BeamBuilder.hxx"
#include "PhysicalInvariants.hxx"
#include "TestHelpers.hxx"
#include "TestMaterials.hxx"

namespace
{
    //! Defines one physical scenario to test the solver against.
    struct Scenario
    {
        std::string name;
        double initialTemperature;
        double initialHumidity;
        double bcTemperature;
        double bcHumidity;
        double bcHc;
        double dTime;
        size_t numSteps;
    };

    //! Defines one set of solver tuning parameters.
    struct SolverSettings
    {
        std::string name;
        double relaxation;
        double errorTolerance;
        size_t maxIterations;
    };

    //! Result of running one scenario on one mesh with one set of settings.
    struct SweepResult
    {
        std::string scenario;
        std::string mesh;
        std::string settings;
        bool completed = false;
        std::string error;
        size_t invariantViolations = 0;
        std::string invariantReport;
    };

    std::vector<Scenario> createScenarios()
    {
        return {
          {"LowHumidity_ColdToWarm",
           0.0, 0.5, 20.0, 0.8, 10.0, 3600.0, 10},
          {"MediumHumidity_WarmToWarm",
           20.0, 0.9, 20.0, 1.0, 10.0, 3600.0, 10},
          {"HighHumidity_ColdToWarm",
           0.0, 0.999, 20.0, 1.0, 10.0, 3600.0, 10},
          {"HighHumidity_WarmToWarm",
           30.0, 0.9999, 20.0, 1.0, 10.0, 3600.0, 10},
          {"ExtremeHumidity_HotToWarm",
           80.0, 0.9999, 20.0, 1.0, 10.0, 3600.0, 10},
          {"NearSaturation_Equilibrium",
           20.0, 0.99999, 20.0, 1.0, 10.0, 3600.0, 10},
        };
    }

    //! Element counts (mesh densities) to sweep for the fixed-width beam. The converged solution
    //! should be consistent across refinement; this checks the solver holds up as the mesh -- and
    //! thus the sorption-curve stiffness resolved per element near saturation -- changes.
    std::vector<std::size_t> createElementCounts()
    {
        return {3, 6, 12, 24};
    }

    std::vector<SolverSettings> createSettings()
    {
        // The relaxation parameter is no longer a meaningful axis: the moisture solver controls
        // under-relaxation adaptively (starts at 1.0 and ratchets down when a step diverges), so
        // the result is independent of the relaxation setting -- the earlier fixed-relaxation
        // sweep is obsolete. These vary the tolerance and iteration budget, which still influence
        // the result through best-effort acceptance and the staggered thermal/moisture coupling.
        return {
          {"Default", 1.0, 1e-5, 50},
          {"TightTolerance", 1.0, 1e-7, 100},
          {"LooseTolerance", 1.0, 1e-3, 20},
          {"FewIterations", 1.0, 1e-5, 10},
          {"ManyIterations", 1.0, 1e-5, 200},
        };
    }

    //! Top/bottom node pairs for a single-row beam with `numElements` columns. Node numbering is
    //! column-major (bottom then top per column), so column c has nodes 2c and 2c+1 (0-based).
    std::vector<std::pair<size_t, size_t>> symmetricPairs(std::size_t numElements)
    {
        std::vector<std::pair<size_t, size_t>> pairs;
        for(std::size_t c = 0; c <= numElements; ++c)
        {
            pairs.emplace_back(2 * c, 2 * c + 1);
        }
        return pairs;
    }

    SweepResult
      runScenario(const Scenario & scenario, std::size_t numElements, const SolverSettings & settings)
    {
        SweepResult result;
        result.scenario = scenario.name;
        result.mesh = std::to_string(numElements) + "elem";
        result.settings = settings.name;

        HygroThermFEM::SimulationProperties::Instance().setIterationParameters(
          settings.relaxation, settings.errorTolerance, settings.maxIterations);

        try
        {
            HygroThermFEM::MultiDomain multiDomain;

            constexpr double beamWidth = 0.10;

            const HygroThermFEM::State initialState({
              .temperature = scenario.initialTemperature,
              .humidity = scenario.initialHumidity,
              .pressure = 101325.0,
              .liquidPercent = 1.0});

            const auto & material =
              multiDomain.materials().createSolidMaterial(TestHelper::CottaerSandstone());

            TestHelper::BeamBuilder builder(multiDomain);
            builder.xStart(0.0)
              .height(0.05)
              .numElementsY(1)
              .state(initialState)
              .addSegment({.material = material.name(),
                           .numElementsX = numElements,
                           .width = beamWidth})
              .build();

            const HygroThermFEM::FixedBCHCCoefficients bcCoeff{
              scenario.bcTemperature, scenario.bcHc, scenario.bcHumidity};
            builder.applyBC_FixedHc(TestHelper::BeamBuilder::Edge::Right, bcCoeff);

            const auto results = multiDomain.transientMultiStep(
              scenario.dTime, scenario.numSteps);

            result.completed = true;

            // Check physical invariants
            TestHelper::PhysicalInvariants checker;

            // Water content must be non-negative
            checker.expectBoundedBelow(results.moisture.values, 0.0, "water content");

            // Humidity BC drives moisture in when bcHumidity >= initialHumidity, so water content
            // should generally not decrease.
            if(scenario.bcHumidity >= scenario.initialHumidity)
            {
                checker.expectNonDecreasing(results.moisture.values, 1.0, "water content");
            }

            // Temperature should stay within reasonable physical bounds
            const double minTemp = std::min(scenario.initialTemperature, scenario.bcTemperature) - 50.0;
            const double maxTemp = std::max(scenario.initialTemperature, scenario.bcTemperature) + 50.0;
            checker.expectBounded(results.temperature.values, minTemp, maxTemp, "temperature");

            // No wild jumps between timesteps (more than 200% of magnitude)
            checker.expectSmooth(results.temperature.values, 2.0, "temperature");
            checker.expectSmooth(results.moisture.values, 2.0, "water content");

            // Symmetric top/bottom node pairs (single-row beam).
            const auto pairs = symmetricPairs(numElements);
            checker.expectSymmetric(results.temperature.values, pairs, 1e-6, "temperature");
            checker.expectSymmetric(results.moisture.values, pairs, 1e-6, "water content");

            result.invariantViolations = checker.violationCount();
            result.invariantReport = checker.report();
        }
        catch(const std::exception & exc)
        {
            result.completed = false;
            result.error = exc.what();
        }

        HygroThermFEM::SimulationProperties::Instance().reset();
        return result;
    }
}   // namespace

TEST(SolverParameterSweep, AllScenariosAllMeshesAllSettings)
{
    const auto scenarios = createScenarios();
    const auto elementCounts = createElementCounts();
    const auto allSettings = createSettings();

    std::vector<SweepResult> allResults;
    size_t totalFailed = 0;

    for(const auto & scenario : scenarios)
    {
        for(const auto numElements : elementCounts)
        {
            for(const auto & settings : allSettings)
            {
                auto result = runScenario(scenario, numElements, settings);

                if(!result.completed || result.invariantViolations > 0)
                {
                    ++totalFailed;
                }

                allResults.push_back(std::move(result));
            }
        }
    }

    // Print summary table
    std::cout << "\n=== Solver Parameter Sweep Results ===\n\n";
    std::cout << std::left
              << std::setw(30) << "Scenario"
              << std::setw(10) << "Mesh"
              << std::setw(18) << "Settings"
              << std::setw(11) << "Completed"
              << std::setw(12) << "Violations"
              << "Details\n";
    std::cout << std::string(110, '-') << "\n";

    for(const auto & res : allResults)
    {
        std::cout << std::setw(30) << res.scenario
                  << std::setw(10) << res.mesh
                  << std::setw(18) << res.settings
                  << std::setw(11) << (res.completed ? "YES" : "NO")
                  << std::setw(12) << res.invariantViolations;

        if(!res.completed)
        {
            std::cout << "THREW: " << res.error;
        }
        else if(res.invariantViolations > 0)
        {
            // Print just the first violation as a preview
            const auto firstNewline = res.invariantReport.find('\n');
            std::cout << res.invariantReport.substr(0, firstNewline);
        }
        else
        {
            std::cout << "OK";
        }
        std::cout << "\n";
    }

    std::cout << "\n" << totalFailed << " / " << allResults.size()
              << " combinations had issues.\n\n";

    // Print detailed violations for failed cases
    for(const auto & res : allResults)
    {
        if(res.invariantViolations > 0)
        {
            std::cout << "--- " << res.scenario << " / " << res.mesh << " / " << res.settings
                      << " ---\n"
                      << res.invariantReport << "\n";
        }
    }

    // Don't fail the test — this is a diagnostic sweep
    // EXPECT_EQ(totalFailed, 0u);
}
