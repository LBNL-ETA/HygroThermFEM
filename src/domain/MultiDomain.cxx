#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <ostream>
#include <ranges>
#include <stdexcept>
#include <utility>

#include "lbnl/algorithm.hxx"

#include "MultiDomain.hxx"
#include "Common.hxx"
#include "FEMMath.hxx"
#include "Functions.hxx"
#include "Nodes.hxx"
#include "MaterialDataChecker.hxx"
#include "SimulationProperties.hxx"

namespace HygroThermFEM
{
    namespace
    {
        //! MultiDomain is the exception boundary for the solver: a single domain reports
        //! non-convergence via ExpectedExt, and here it is turned back into the runtime_error the
        //! application layer already expects and reports to the user.
        SingleSolution solveOrThrow(const lbnl::ExpectedExt<SingleSolution, SolverError> & result)
        {
            if(!result.has_value())
            {
                throw std::runtime_error("Solution failed to converge.");
            }
            return result.value();
        }

        //! Liquid fractions lambda(T) per node for the freezing bookkeeping.
        std::vector<double> liquidPercentsFromTemperatures(const std::vector<double> & temps)
        {
            return lbnl::transform_to_vector(
              temps, [](const double temp) { return liquidFractionFromTemperature(temp); });
        }
    }   // namespace

    MultiDomain::MultiDomain() : MultiDomain(MultiDomainParams{})
    {}

    // passing false to subdomains means that previous timestep values will not be automatically
    // updated. This mean that multidomain must update its values once solution converged.
    MultiDomain::MultiDomain(MultiDomainParams params) :
        m_ThermalDomain(m_Nodes, m_Materials, /*automaticUpdatePreviousTimestep=*/false),
        m_MoistureDomain(m_Nodes, m_Materials, /*automaticUpdatePreviousTimestep=*/false),
        m_SimulateThermal(params.performThermal),
        m_SimulateMoisture(params.performMoisture)
    {}

    Solution MultiDomain::transient(const std::vector<double> & temperature,
                                    const std::vector<double> & humidity,
                                    const double t_DTime,
                                    const size_t timestepIndex)
    {
        auto currentTemperature{temperature};
        auto currentHumidity{humidity};
        double totalTime{0};
        size_t subStep{0};

        auto temperatureError{std::numeric_limits<double>::max()};
        auto humidityError{std::numeric_limits<double>::max()};

        // Adaptive time-stepping loop: accumulate sub-steps until the full
        // requested dt is covered.  Each domain returns after its first
        // successful solve (possibly at a smaller dt), and cross-coupling
        // data is exchanged after every sub-step.  This avoids grinding
        // thousands of internal sub-steps with frozen cross-domain fields.
        while(totalTime < t_DTime)
        {
            const double remainingTime = t_DTime - totalTime;
            double effectiveDt = remainingTime;

            if(m_DiagStream != nullptr)
            {
                *m_DiagStream << "\n# substep=" << subStep
                              << " totalTime=" << std::scientific << std::setprecision(6)
                              << totalTime
                              << " remaining=" << remainingTime << "\n";
            }

            // Solve moisture first (typically the stiff domain).
            // It may return a smaller dt than requested.
            SingleSolution humiditySolution{currentHumidity, effectiveDt};
            if(m_SimulateMoisture)
            {
                if(m_DiagStream != nullptr)
                {
                    *m_DiagStream << "\n# MOISTURE substep=" << subStep << "\n";
                }
                humiditySolution =
                  solveOrThrow(m_MoistureDomain.transient(currentHumidity, effectiveDt, timestepIndex));
                effectiveDt = std::min(effectiveDt, humiditySolution.dTime);
            }

            // Solve thermal at the (possibly reduced) dt.
            SingleSolution temperatureSolution{currentTemperature, effectiveDt};
            if(m_SimulateThermal)
            {
                if(m_DiagStream != nullptr)
                {
                    *m_DiagStream << "\n# THERMAL substep=" << subStep
                                  << " dt=" << effectiveDt << "\n";
                }
                temperatureSolution =
                  solveOrThrow(m_ThermalDomain.transient(currentTemperature, effectiveDt, timestepIndex));

                // If thermal needed an even smaller dt, redo moisture at that dt.
                if(temperatureSolution.dTime < effectiveDt && m_SimulateMoisture)
                {
                    effectiveDt = temperatureSolution.dTime;
                    if(m_DiagStream != nullptr)
                    {
                        *m_DiagStream << "\n# MOISTURE redo at dt=" << effectiveDt << "\n";
                    }
                    humiditySolution =
                      solveOrThrow(m_MoistureDomain.transient(currentHumidity, effectiveDt, timestepIndex));
                }
            }

            // Compute convergence errors for reporting
            humidityError = normError(humiditySolution.solution, currentHumidity);
            temperatureError = normError(temperatureSolution.solution, currentTemperature);

            // Accept solutions and exchange cross-coupling data.
            currentHumidity = humiditySolution.solution;
            currentTemperature = temperatureSolution.solution;
            m_Nodes.updateNodeTemperatures(currentTemperature, false);
            m_Nodes.updateNodeHumidities(currentHumidity, false);
            if(!SimulationProperties::Instance().excludeLatentHeatOfFusion())
            {
                // Roll the freezing bookkeeping with the accepted temperatures so the
                // liquid/ice split the capacitance terms read stays consistent.
                m_Nodes.updateNodeLiquidPercents(
                  liquidPercentsFromTemperatures(currentTemperature), false);
            }

            totalTime += effectiveDt;
            ++subStep;

            if(totalTime < t_DTime)
            {
                // More substeps follow: advance the nodes' previous-timestep values to the
                // accepted substep state so the secant storage capacities (sorption,
                // fusion) span exactly the interval the mass matrices book on the next
                // substep, [substep start, iterate]. Leaving previous at the full step's
                // start would leak storage energy at every substep boundary (same
                // telescoping argument as IDomain::resolveShockStep). The current values
                // were just set above, so a single previous-updating call rolls them.
                m_Nodes.updateNodeTemperatures(currentTemperature, true);
                m_Nodes.updateNodeHumidities(currentHumidity, true);
                if(!SimulationProperties::Instance().excludeLatentHeatOfFusion())
                {
                    m_Nodes.updateNodeLiquidPercents(
                      liquidPercentsFromTemperatures(currentTemperature), true);
                }
            }

            if(m_DiagStream != nullptr)
            {
                *m_DiagStream << "\n# accepted dt=" << std::scientific << std::setprecision(6)
                              << effectiveDt
                              << " totalTime=" << totalTime
                              << " temp_err=" << temperatureError
                              << " hum_err=" << humidityError << "\n";
            }
        }

        // Coupling convergence check: with the latest cross-coupling data
        // baked into nodes, re-solve each domain for one tiny step from the
        // current state. If the coupled solution is self-consistent, the
        // solution should barely change. The reported error measures the
        // coupling residual — how much the solution shifts when both domains
        // see each other's latest fields.
        if(m_DiagStream != nullptr)
        {
            *m_DiagStream << "\n# CONVERGENCE CHECK\n";
        }

        // Use a small probe dt — just enough for the NR to converge and
        // reveal any coupling mismatch, without advancing significantly.
        const double probeDt = t_DTime * 1e-6;

        if(m_SimulateMoisture)
        {
            if(m_DiagStream != nullptr)
            {
                *m_DiagStream << "\n# MOISTURE convergence probe\n";
            }
            m_Nodes.updateNodeTemperatures(currentTemperature, false);
            const auto humProbe =
              solveOrThrow(m_MoistureDomain.transient(currentHumidity, probeDt, timestepIndex));
            humidityError = normError(humProbe.solution, currentHumidity);
        }

        if(m_SimulateThermal)
        {
            if(m_DiagStream != nullptr)
            {
                *m_DiagStream << "\n# THERMAL convergence probe\n";
            }
            m_Nodes.updateNodeHumidities(currentHumidity, false);
            const auto tempProbe =
              solveOrThrow(m_ThermalDomain.transient(currentTemperature, probeDt, timestepIndex));
            temperatureError = normError(tempProbe.solution, currentTemperature);
        }

        if(m_DiagStream != nullptr)
        {
            *m_DiagStream << "\n# FINAL temp_err=" << std::scientific << std::setprecision(6)
                          << temperatureError
                          << " hum_err=" << humidityError << "\n";
        }

        // Final update: advance the "previous timestep" values in nodes
        m_Nodes.updateNodeTemperatures(currentTemperature, true);
        m_Nodes.updateNodeHumidities(currentHumidity, true);
        if(!SimulationProperties::Instance().excludeLatentHeatOfFusion())
        {
            m_Nodes.updateNodeLiquidPercents(
              liquidPercentsFromTemperatures(currentTemperature), true);
        }

        const auto waterContent = m_Nodes.properties(Variable::water);
        const auto liquidContent = m_Nodes.properties(Variable::liquid);
        const auto vaporContent = m_Nodes.properties(Variable::vapor);
        const auto iceContent = m_Nodes.properties(Variable::ice);

        const auto heatFlux = m_ThermalDomain.flux();
        const auto waterFlux = m_MoistureDomain.flux();

        return Solution{t_DTime,
                        currentTemperature,
                        currentHumidity,
                        waterContent,
                        liquidContent,
                        vaporContent,
                        iceContent,
                        heatFlux,
                        waterFlux,
                        temperatureError,
                        humidityError};
    }

    TransientResults MultiDomain::transientMultiStep(const double dTime, const size_t numSteps)
    {
        auto temperatures = m_Nodes.properties(Variable::temperature);
        auto humidities = m_Nodes.properties(Variable::humidity);

        TransientResults results;

        for(size_t step = 0; step < numSteps; ++step)
        {
            const auto solution = transient(temperatures, humidities, dTime, step);

            results.temperature.values.push_back(solution.temperature);
            results.temperature.errors.push_back(solution.temperatureError);
            results.moisture.values.push_back(solution.waterContent);
            results.moisture.errors.push_back(solution.humidityError);

            temperatures = solution.temperature;
            humidities = solution.humidity;
        }

        return results;
    }

    void MultiDomain::setSolverSettings(const SolverSettings & settings)
    {
        m_SolverSettings = settings;
        m_ThermalDomain.setSolverSettings(settings);
        m_MoistureDomain.setSolverSettings(settings);
    }

    SolverSettings MultiDomain::solverSettings() const
    {
        return m_SolverSettings ? *m_SolverSettings : SolverSettings::fromGlobals();
    }

    Solution MultiDomain::steadyState()
    {
        const auto settings = solverSettings();
        const auto ConvergenceError = settings.errorTolerance;
        auto temperatureError{std::numeric_limits<double>::max()};
        auto humidityError{std::numeric_limits<double>::max()};
        const auto MaxIterations = settings.maxNumberOfIterations;
        size_t currentIteration{0};
        auto humidity = m_Nodes.properties(Variable::humidity);
        auto previousHumidity = humidity;
        auto temperature = m_Nodes.properties(Variable::temperature);
        auto previousTemperature = temperature;
        // Steady-state coupling loop: alternately solve moisture and thermal domains,
        // exchanging cross-domain data after each solve, until both domains converge.
        // Uses || (either unconverged) so the loop only exits when both are satisfied.
        do
        {
            if(m_SimulateMoisture)
            {
                // steadyState() enforces the domain's solutionBounds: condensation zones pin at
                // saturation instead of reporting Glaser-style supersaturation.
                humidity = m_MoistureDomain.steadyState();
                humidityError = normError(humidity, previousHumidity);
                previousHumidity = humidity;
                m_Nodes.updateNodeHumidities(humidity);
            }
            else
            {
                humidityError = 0;
            }
            if(m_SimulateThermal)
            {
                temperature = m_ThermalDomain.steadyState();
                temperatureError = normError(temperature, previousTemperature);
                previousTemperature = temperature;
                m_Nodes.updateNodeTemperatures(temperature);
            }
            else
            {
                temperatureError = 0;
            }
            ++currentIteration;
        } while((temperatureError > ConvergenceError || humidityError > ConvergenceError)
                && currentIteration < MaxIterations);

        m_Nodes.updateNodeHumidities(humidity, true);
        m_Nodes.updateNodeTemperatures(temperature, true);

        const auto waterContent = m_Nodes.properties(Variable::water);
        const auto liquidContent = m_Nodes.properties(Variable::liquid);
        const auto vaporContent = m_Nodes.properties(Variable::vapor);
        const auto iceContent = m_Nodes.properties(Variable::ice);

        // Only compute a domain's flux when that domain is actually simulated. A thermal-only
        // steady run has no moisture material properties, so computing the moisture flux would
        // dereference missing data.
        const auto heatFlux = m_SimulateThermal ? m_ThermalDomain.flux() : std::vector<NodeFlux>{};
        const auto waterFlux = m_SimulateMoisture ? m_MoistureDomain.flux() : std::vector<NodeFlux>{};

        return Solution{0,
                        temperature,
                        humidity,
                        waterContent,
                        liquidContent,
                        vaporContent,
                        iceContent,
                        heatFlux,
                        waterFlux,
                        temperatureError,
                        humidityError};
    }

    Solution MultiDomain::currentStateSolution()
    {
        const auto temperature = m_Nodes.properties(Variable::temperature);
        const auto humidity = m_Nodes.properties(Variable::humidity);
        const auto waterContent = m_Nodes.properties(Variable::water);
        const auto liquidContent = m_Nodes.properties(Variable::liquid);
        const auto vaporContent = m_Nodes.properties(Variable::vapor);
        const auto iceContent = m_Nodes.properties(Variable::ice);

        // No step has been taken from this state, so there is nothing to differentiate: the
        // fluxes are identically zero (exact for the uniform initial condition).
        const std::vector<NodeFlux> zeroFlux(temperature.size(), NodeFlux{0.0, 0.0});

        return Solution{0.0,
                        temperature,
                        humidity,
                        waterContent,
                        liquidContent,
                        vaporContent,
                        iceContent,
                        zeroFlux,
                        zeroFlux,
                        0.0,
                        0.0};
    }

    void MultiDomain::setDiagnosticStream(std::ostream * stream)
    {
        m_DiagStream = stream;
        m_ThermalDomain.setDiagnosticStream(stream);
        m_MoistureDomain.setDiagnosticStream(stream);
    }

    void MultiDomain::performMoistureSimulation(const bool val)
    {
        m_SimulateMoisture = val;
    }

    void MultiDomain::performThermalSimulation(const bool val)
    {
        m_SimulateThermal = val;
    }

    void MultiDomain::createElement(const ElementParams & params)
    {
        // At least one domain must create elements for material assignment to nodes
        // (needed for water content calculations even without active simulations)
        if(m_SimulateThermal || !m_SimulateMoisture)
        {
            m_ThermalDomain.createElement(params.node1, params.node2, params.node3, params.node4, params.material);
        }
        if(m_SimulateMoisture)
        {
            m_MoistureDomain.createElement(params.node1, params.node2, params.node3, params.node4, params.material);
        }
    }

    void MultiDomain::createBC_FixedHc(const size_t index1,
                                       const size_t index2,
                                       const FixedBCHCCoefficients & fixedBchcCoefficients)
    {
        m_ThermalDomain.createBC_FixedHc(index1, index2, fixedBchcCoefficients, m_SimulateMoisture);
        if(m_SimulateMoisture)
        {
            m_MoistureDomain.createBC_FixedHc(index1, index2, fixedBchcCoefficients);
        }
    }

    void MultiDomain::createBC_FixedHc(
      size_t index1,
      size_t index2,
      const std::vector<FixedBCHCCoefficients> & fixedBchcCoefficients)
    {
        m_ThermalDomain.createBC_FixedHc(index1, index2, fixedBchcCoefficients, m_SimulateMoisture);
        if(m_SimulateMoisture)
        {
            m_MoistureDomain.createBC_FixedHc(index1, index2, fixedBchcCoefficients);
        }
    }

    void MultiDomain::createBC_TARPHc(size_t index1,
                                      size_t index2,
                                      const TARPCoefficients & varHCCoeff,
                                      double surfaceTilt)
    {
        m_ThermalDomain.createBC_TARPHc(
          index1, index2, varHCCoeff, surfaceTilt, m_SimulateMoisture);
        if(m_SimulateMoisture)
        {
            m_MoistureDomain.createBC_TARPHc(index1, index2, varHCCoeff, surfaceTilt);
        }
    }

    void MultiDomain::createBC_TARPHc(size_t index1,
                                      size_t index2,
                                      const std::vector<TARPCoefficients> & varHCCoeff,
                                      double surfaceTilt)
    {
        m_ThermalDomain.createBC_TARPHc(
          index1, index2, varHCCoeff, surfaceTilt, m_SimulateMoisture);
        if(m_SimulateMoisture)
        {
            m_MoistureDomain.createBC_TARPHc(index1, index2, varHCCoeff, surfaceTilt);
        }
    }

    void MultiDomain::createBC_ASHRAEInsideHc(size_t index1,
                                              size_t index2,
                                              const ASHRAEInsideCoefficients & coeff,
                                              double surfaceHeight,
                                              double surfaceTilt)
    {
        m_ThermalDomain.createBC_ASHRAEInsideHc(
          index1, index2, coeff, surfaceHeight, surfaceTilt, m_SimulateMoisture);
        if(m_SimulateMoisture)
        {
            m_MoistureDomain.createBC_ASHRAEInsideHc(index1, index2, coeff, surfaceHeight, surfaceTilt);
        }
    }

    void MultiDomain::createBC_ASHRAEInsideHc(size_t index1,
                                              size_t index2,
                                              const std::vector<ASHRAEInsideCoefficients> & coeff,
                                              double surfaceHeight,
                                              double surfaceTilt)
    {
        m_ThermalDomain.createBC_ASHRAEInsideHc(
          index1, index2, coeff, surfaceHeight, surfaceTilt, m_SimulateMoisture);
        if(m_SimulateMoisture)
        {
            m_MoistureDomain.createBC_ASHRAEInsideHc(index1, index2, coeff, surfaceHeight, surfaceTilt);
        }
    }

    void MultiDomain::createBC_ASHRAEOutsideHc(size_t index1,
                                               size_t index2,
                                               const ASHRAEOutsideCoefficients & coeff)
    {
        m_ThermalDomain.createBC_ASHRAEOutsideHc(index1, index2, coeff, m_SimulateMoisture);
        if(m_SimulateMoisture)
        {
            m_MoistureDomain.createBC_ASHRAEOutsideHc(index1, index2, coeff);
        }
    }

    void MultiDomain::createBC_ASHRAEOutsideHc(size_t index1,
                                               size_t index2,
                                               const std::vector<ASHRAEOutsideCoefficients> & coeff)
    {
        m_ThermalDomain.createBC_ASHRAEOutsideHc(index1, index2, coeff, m_SimulateMoisture);
        if(m_SimulateMoisture)
        {
            m_MoistureDomain.createBC_ASHRAEOutsideHc(index1, index2, coeff);
        }
    }

    void MultiDomain::createBC_YazdanianKlemsHc(size_t index1,
                                                size_t index2,
                                                const YazdanianKlemsCoefficients & coeff)
    {
        m_ThermalDomain.createBC_YazdanianKlemsHc(index1, index2, coeff, m_SimulateMoisture);
        if(m_SimulateMoisture)
        {
            m_MoistureDomain.createBC_YazdanianKlemsHc(index1, index2, coeff);
        }
    }

    void MultiDomain::createBC_YazdanianKlemsHc(
      size_t index1, size_t index2, const std::vector<YazdanianKlemsCoefficients> & coeff)
    {
        m_ThermalDomain.createBC_YazdanianKlemsHc(index1, index2, coeff, m_SimulateMoisture);
        if(m_SimulateMoisture)
        {
            m_MoistureDomain.createBC_YazdanianKlemsHc(index1, index2, coeff);
        }
    }

    void
      MultiDomain::createBC_KimuraHc(size_t index1, size_t index2, const KimuraCoefficients & coeff)
    {
        m_ThermalDomain.createBC_KimuraHc(index1, index2, coeff, m_SimulateMoisture);
        if(m_SimulateMoisture)
        {
            m_MoistureDomain.createBC_KimuraHc(index1, index2, coeff);
        }
    }

    void MultiDomain::createBC_KimuraHc(size_t index1,
                                        size_t index2,
                                        const std::vector<KimuraCoefficients> & coeff)
    {
        m_ThermalDomain.createBC_KimuraHc(index1, index2, coeff, m_SimulateMoisture);
        if(m_SimulateMoisture)
        {
            m_MoistureDomain.createBC_KimuraHc(index1, index2, coeff);
        }
    }

    void MultiDomain::createBC_FixedTemperature(const size_t index1,
                                                const size_t index2,
                                                const double t_Temp1,
                                                const double t_Temp2)
    {
        m_ThermalDomain.createBC_FixedTemperature(index1, index2, t_Temp1, t_Temp2);
    }

    void MultiDomain::createBC_FixedTemperature(size_t index1,
                                                size_t index2,
                                                const std::vector<ConstantBCTemperatures> & temp)
    {
        m_ThermalDomain.createBC_FixedTemperature(index1, index2, temp);
    }

    void MultiDomain::createBC_FixedTemperature(const size_t index1,
                                                const size_t index2,
                                                const double t_Temp)
    {
        m_ThermalDomain.createBC_FixedTemperature(index1, index2, t_Temp);
    }

    void
      MultiDomain::createBC_FixedTemperature(size_t index1, size_t index2, std::vector<double> temp)
    {
        m_ThermalDomain.createBC_FixedTemperature(index1, index2, std::move(temp));
    }

    void MultiDomain::createBC_FixedTemperatureAndHumidity(size_t index1,
                                                           size_t index2,
                                                           const TemperatureAndHumidity & values)
    {
        m_ThermalDomain.createBC_FixedTemperature(index1, index2, values.Temperature);
        if(m_SimulateMoisture)
        {
            m_MoistureDomain.createBC_FixedHumidity(index1, index2, values);
        }
    }

    void MultiDomain::createBC_FixedTemperatureAndHumidity(
      size_t index1, size_t index2, const std::vector<TemperatureAndHumidity> & values)
    {
        std::vector<double> temperatures;
        temperatures.reserve(values.size());
        std::ranges::transform(values, std::back_inserter(temperatures),
            [](const TemperatureAndHumidity & val) { return val.Temperature; });
        m_ThermalDomain.createBC_FixedTemperature(index1, index2, temperatures);
        if(m_SimulateMoisture)
        {
            m_MoistureDomain.createBC_FixedHumidity(index1, index2, values);
        }
    }

    void MultiDomain::createBC_FixedHeatFlux(size_t index1, size_t index2, double t_Flux)
    {
        m_ThermalDomain.createBC_FixedFlux(index1, index2, t_Flux);
    }

    void
      MultiDomain::createBC_FixedHeatFlux(size_t index1, size_t index2, std::vector<double> t_Flux)
    {
        m_ThermalDomain.createBC_FixedFlux(index1, index2, t_Flux);
    }

    void MultiDomain::createBC_BlackBodyRadiation(size_t index1,
                                                  size_t index2,
                                                  double t_Emissivity,
                                                  double t_RadiationTemperature)
    {
        m_ThermalDomain.createBC_BlackBodyRadiation(
          index1, index2, t_Emissivity, t_RadiationTemperature);
    }

    void MultiDomain::createBC_BlackBodyRadiation(
      size_t index1, size_t index2, const std::vector<BlackBodyRadiationBCCoefficients> & radCoeffs)
    {
        m_ThermalDomain.createBC_BlackBodyRadiation(index1, index2, radCoeffs);
    }

    void MultiDomain::createBC_LinearizedRadiation(
      const size_t index1,
      const size_t index2,
      const LinearizedRadiationBCCoefficients & linearRadBC)
    {
        m_ThermalDomain.createBC_LinearizedRadiation(index1, index2, linearRadBC);
    }

    void MultiDomain::createBC_LinearizedRadiation(
      size_t index1,
      size_t index2,
      const std::vector<LinearizedRadiationBCCoefficients> & linearRadBC)
    {
        m_ThermalDomain.createBC_LinearizedRadiation(index1, index2, linearRadBC);
    }

    void MultiDomain::createEnclosureRadiation(
      const std::vector<EnclosureRadiationSegment> & segments,
      const std::map<std::size_t, double> & openEnclosureTemperatures,
      const bool smoothViewFactors,
      const EnclosureSurfaceTemperature surfaceTemperature)
    {
        m_ThermalDomain.createEnclosureRadiation(
          segments, openEnclosureTemperatures, smoothViewFactors, surfaceTemperature);
    }

    double MultiDomain::normError(const std::vector<double> & vec1,
                                  const std::vector<double> & vec2)
    {
        auto norm1 = norm(vec1);
        auto norm2 = norm(vec2);
        if(norm1 == 0)
        {
            norm1 = 1e-10;
            if(norm2 == 0)
            {
                norm2 = norm1;
            }
        }

        return std::abs(norm1 - norm2) / norm1;
    }

    std::vector<double> MultiDomain::property(const Variable property) const
    {
        return m_Nodes.properties(property);
    }

    void MultiDomain::setGravityVector(const FenestrationCommon::GravityVector & gravityVector)
    {
        m_ThermalDomain.setGravityVector(gravityVector);
    }

    void MultiDomain::subscribeThermal(Timesteps::TimestepObserver * observer)
    {
        m_ThermalDomain.subscribe(observer);
    }

    void MultiDomain::unsubscribeThermal(Timesteps::TimestepObserver * observer)
    {
        m_ThermalDomain.unsubscribe(observer);
    }

    void MultiDomain::subscribeMoisture(Timesteps::TimestepObserver * observer)
    {
        m_MoistureDomain.subscribe(observer);
    }

    void MultiDomain::unsubscribeMoisture(Timesteps::TimestepObserver * observer)
    {
        m_MoistureDomain.unsubscribe(observer);
    }

    bool MultiDomain::isMoistureSimulationON() const
    {
        return m_SimulateMoisture;
    }

    bool MultiDomain::isThermalSimulationON() const
    {
        return m_SimulateThermal;
    }

    MaterialsErrorCheckVector MultiDomain::checkMaterialsForTransientSimulation() const
    {
        MaterialDataChecker dataChecker{*this};
        return dataChecker.checkMaterialProperties(true);
    }

    MaterialsErrorCheckVector MultiDomain::checkMaterialsForSteadyStateSimulation() const
    {
        // The full flag-aware check: with moisture off it reduces to the thermal essentials
        // (dry conductivity + emissivity), and with moisture on it requires the moisture
        // properties too -- a steady moisture run on hole-riddled legacy materials must be
        // stopped with a report instead of reaching the engine.
        MaterialDataChecker dataChecker{*this};
        return dataChecker.checkMaterialProperties(false);
    }

    MaterialsErrorCheckVector
      MultiDomain::checkForMaterialsValidity(const SimulationType simulationType) const
    {
        MaterialsErrorCheckVector result;
        switch(simulationType)
        {
            case SimulationType::SteadyState:
                result = checkMaterialsForSteadyStateSimulation();
                break;
            case SimulationType::Transient:
                result = checkMaterialsForTransientSimulation();
                break;
            default:
                throw std::runtime_error("Incorrect selection of simulation type.");
        }
        return result;
    }

    void MultiDomain::clearModel()
    {
        m_ThermalDomain.clearModel();
        m_MoistureDomain.clearModel();
    }

    Materials & MultiDomain::materials()
    {
        return m_Materials;
    }

    const Materials & MultiDomain::materials() const
    {
        return m_Materials;
    }

    ThermalDomain & MultiDomain::thermal()
    {
        return m_ThermalDomain;
    }

    MoistureDomain & MultiDomain::moisture()
    {
        return m_MoistureDomain;
    }

    Nodes & MultiDomain::nodes()
    {
        return m_Nodes;
    }

    const Nodes & MultiDomain::nodes() const
    {
        return m_Nodes;
    }

    Solution::Solution(const double dtime,
                       std::vector<double> temperature,
                       std::vector<double> humidity,
                       std::vector<double> waterContent,
                       std::vector<double> liquidWaterContent,
                       std::vector<double> vaporContent,
                       std::vector<double> iceContent,
                       std::vector<NodeFlux> heatFlux,
                       std::vector<NodeFlux> waterFlux,
                       const double temperatureError,
                       const double humidityError) :
        dTime(dtime),
        temperature(std::move(temperature)),
        humidity(std::move(humidity)),
        waterContent(std::move(waterContent)),
        liquidWaterContent(std::move(liquidWaterContent)),
        vaporContent(std::move(vaporContent)),
        iceContent(std::move(iceContent)),
        heatFlux(std::move(heatFlux)),
        waterFlux(std::move(waterFlux)),
        temperatureError(temperatureError),
        humidityError(humidityError)
    {}
}   // namespace HygroThermFEM
