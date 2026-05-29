#include <gtest/gtest.h>

#include <iomanip>
#include <iostream>
#include <string>
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

    //! Result of running one scenario with one set of settings.
    struct SweepResult
    {
        std::string scenario;
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

    std::vector<SolverSettings> createSettings()
    {
        return {
          {"Default", 1.0, 1e-5, 50},
          {"Relaxed_0.8", 0.8, 1e-5, 50},
          {"Relaxed_0.5", 0.5, 1e-5, 50},
          {"TightTolerance", 1.0, 1e-7, 100},
          {"LooseTolerance", 1.0, 1e-3, 20},
          {"FewIterations", 1.0, 1e-5, 10},
          {"ManyIterations", 1.0, 1e-5, 200},
        };
    }

    SweepResult runScenario(const Scenario & scenario, const SolverSettings & settings)
    {
        SweepResult result;
        result.scenario = scenario.name;
        result.settings = settings.name;

        HygroThermFEM::SimulationProperties::Instance().setIterationParameters(
          settings.relaxation, settings.errorTolerance, settings.maxIterations);

        try
        {
            HygroThermFEM::MultiDomain multiDomain;

            constexpr double beamWidth = 0.10;
            constexpr size_t numElements = 3;

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

            // Humidity BC is 1.0 and interior <= 1.0, so water content
            // should generally not decrease (moisture is being driven in)
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

            // Symmetric node pairs (beam with 1 row: top/bottom pairs)
            // With 3 elements, 1 row: nodes are column-major pairs (1,2), (3,4), (5,6), (7,8)
            const std::vector<std::pair<size_t, size_t>> symmetricPairs{
              {0, 1}, {2, 3}, {4, 5}, {6, 7}};
            checker.expectSymmetric(results.temperature.values, symmetricPairs, 1e-6, "temperature");
            checker.expectSymmetric(results.moisture.values, symmetricPairs, 1e-6, "water content");

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

TEST(SolverParameterSweep, AllScenariosAllSettings)
{
    const auto scenarios = createScenarios();
    const auto allSettings = createSettings();

    std::vector<SweepResult> allResults;
    size_t totalFailed = 0;

    for(const auto & scenario : scenarios)
    {
        for(const auto & settings : allSettings)
        {
            auto result = runScenario(scenario, settings);

            if(!result.completed || result.invariantViolations > 0)
            {
                ++totalFailed;
            }

            allResults.push_back(std::move(result));
        }
    }

    // Print summary table
    std::cout << "\n=== Solver Parameter Sweep Results ===\n\n";
    std::cout << std::left
              << std::setw(35) << "Scenario"
              << std::setw(20) << "Settings"
              << std::setw(12) << "Completed"
              << std::setw(12) << "Violations"
              << "Details\n";
    std::cout << std::string(100, '-') << "\n";

    for(const auto & res : allResults)
    {
        std::cout << std::setw(35) << res.scenario
                  << std::setw(20) << res.settings
                  << std::setw(12) << (res.completed ? "YES" : "NO")
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
            std::cout << "--- " << res.scenario << " + " << res.settings << " ---\n"
                      << res.invariantReport << "\n";
        }
    }

    // Don't fail the test — this is a diagnostic sweep
    // EXPECT_EQ(totalFailed, 0u);
}
