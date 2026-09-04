#pragma once

#include <cmath>
#include <cstddef>
#include <vector>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

//! Shared machinery for the SteadyVsTransient_* tests. Each of those files supplies a wall;
//! the method below is the same for all of them.
//!
//! WHAT THESE TESTS ARE FOR. Under constant boundary conditions a transient run that has
//! stopped changing IS a steady solution of the same equations. These tests assert exactly
//! that, for every moisture modelling option in turn, and report how far apart the two
//! solvers have drifted when they fail.
//!
//! Nothing is stored. The transient endpoint is computed in the test and the steady solve is
//! compared against it directly, so the check depends on no recorded number, needs no
//! re-baselining, and keeps its meaning after any deliberate change to the engine. Pinning
//! absolute values is the job of the golden tests elsewhere in this suite.
//!
//! The cost is one march to rest per case, a few seconds. How long that march is, and the
//! separation the two solvers are allowed, are both properties of the wall, so each test file
//! states its own schedule and tolerances rather than inheriting hidden ones from here.
namespace TestHelper::SteadyVsTransient
{
    //! The march used to bring the transient to rest. Each test file states its own, so a
    //! reader can see what was actually run without going looking for it.
    struct TransientSchedule
    {
        std::size_t steps{0u};
        double stepSeconds{0.0};

        [[nodiscard]] double simulatedYears() const
        {
            constexpr double secondsPerYear{365.0 * 24.0 * 3600.0};
            return static_cast<double>(steps) * stepSeconds / secondsPerYear;
        }
    };

    //! Fields per x-column, exterior face first. Walls here are one element tall, so the two
    //! nodes of a column carry the same value.
    struct Fields
    {
        std::vector<double> temperature;   //!< C
        std::vector<double> humidity;      //!< RH fraction
    };

    //! The moisture modelling options, named and defaulted the way the THERM dialog presents
    //! them: every box starts unticked, and a case lists only what it switches on.
    struct Options
    {
        bool liquidTransport{false};
        bool heatOfEvaporation{false};
        bool capillaryConduction{false};
        bool vaporDiffusionConduction{false};
    };

    //! How far the two solvers may sit apart. This is not solver error but the transient's
    //! remaining distance from rest, so it belongs to the case: how fast a wall settles
    //! depends on which options are switched on, and the cases below differ by orders of
    //! magnitude. Each states its own measured value.
    struct Tolerances
    {
        double temperature{1e-3};   //!< C
        double humidity{1e-3};      //!< RH fraction
    };

    struct ModellingCase
    {
        const char * name{""};
        Options options{};
        //! Ten simulated years, the default march. Lengthening it is the lever when a case
        //! needs a tighter tolerance than it currently earns.
        TransientSchedule schedule{.steps = 8760u, .stepSeconds = 10.0 * 3600.0};
        Tolerances tolerances{};
    };

    //! The engine takes these as EXCLUDE flags, so the inversion happens here, once.
    inline void applyOptions(const Options & options)
    {
        HygroThermFEM::SimulationProperties::Instance().setCalculationParameters(
          !options.liquidTransport,
          !options.heatOfEvaporation,
          !options.capillaryConduction,
          !options.vaporDiffusionConduction,
          false);
    }

    inline std::vector<double> perColumn(const std::vector<double> & nodal)
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
        Fields fields;
        bool converged{false};
        std::size_t passes{0u};
    };

    template<typename BuildWall>
    SteadyRun runSteady(const BuildWall & build, const ModellingCase & modellingCase)
    {
        applyOptions(modellingCase.options);
        HygroThermFEM::MultiDomain multiDomain;
        build(multiDomain);
        const auto solution = multiDomain.steadyState();
        HygroThermFEM::SimulationProperties::Instance().resetCalculationParameters();
        return {{perColumn(solution.temperature), perColumn(solution.humidity)},
                solution.converged,
                solution.steadyPasses};
    }

    template<typename BuildWall>
    Fields runTransient(const BuildWall & build, const ModellingCase & modellingCase)
    {
        const auto & schedule = modellingCase.schedule;
        applyOptions(modellingCase.options);
        HygroThermFEM::MultiDomain multiDomain;
        build(multiDomain);
        auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
        auto humidities = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);
        for(std::size_t step = 0u; step < schedule.steps; ++step)
        {
            const auto solution =
              multiDomain.transient(temperatures, humidities, schedule.stepSeconds, step);
            temperatures = solution.temperature;
            humidities = solution.humidity;
        }
        HygroThermFEM::SimulationProperties::Instance().resetCalculationParameters();
        return {perColumn(temperatures), perColumn(humidities)};
    }

    struct Drift
    {
        double worst{0.0};
        std::size_t column{0u};
    };

    inline Drift worstDrift(const std::vector<double> & steady,
                            const std::vector<double> & transient)
    {
        Drift drift;
        for(std::size_t column = 0u; column < transient.size(); ++column)
        {
            const double separation{std::abs(steady[column] - transient[column])};
            if(separation > drift.worst)
            {
                drift.worst = separation;
                drift.column = column;
            }
        }
        return drift;
    }

    //! One assertion per field, on the worst separation anywhere in the wall. A failure then
    //! leads with the number the reader wants -- how far apart the two solvers have drifted
    //! -- and says where, in what, and over how long a run. A wall that agrees says nothing
    //! at all, which is why the drift is not printed: nobody reads a passing test.
    inline void expectDriftWithin(const std::vector<double> & steady,
                                  const std::vector<double> & transient,
                                  const double tolerance,
                                  const char * caseName,
                                  const char * quantity,
                                  const TransientSchedule & schedule)
    {
        ASSERT_EQ(steady.size(), transient.size()) << caseName << ", " << quantity;
        const auto drift = worstDrift(steady, transient);
        EXPECT_LE(drift.worst, tolerance)
          << caseName << ": steady and transient " << quantity << " are " << drift.worst
          << " apart at column " << drift.column + 1u << " of " << transient.size() << " (steady "
          << steady[drift.column] << ", transient " << transient[drift.column] << "), after "
          << schedule.steps << " steps of " << schedule.stepSeconds / 3600.0 << " h, "
          << schedule.simulatedYears() << " simulated years";
    }

    //! The whole test: march the transient to rest, solve the same case steadily, and require
    //! them to be the same field.
    template<typename BuildWall>
    void expectSteadyMatchesTransient(const BuildWall & build, const ModellingCase & modellingCase)
    {
        const auto transient = runTransient(build, modellingCase);
        const auto steady = runSteady(build, modellingCase);

        EXPECT_TRUE(steady.converged)
          << modellingCase.name << ": steady solve stopped on its pass budget after "
          << steady.passes << " passes";
        expectDriftWithin(steady.fields.temperature,
                          transient.temperature,
                          modellingCase.tolerances.temperature,
                          modellingCase.name,
                          "temperature",
                          modellingCase.schedule);
        expectDriftWithin(steady.fields.humidity,
                          transient.humidity,
                          modellingCase.tolerances.humidity,
                          modellingCase.name,
                          "humidity",
                          modellingCase.schedule);
    }
}   // namespace TestHelper::SteadyVsTransient
