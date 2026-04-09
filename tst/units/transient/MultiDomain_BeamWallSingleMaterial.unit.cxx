#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

#include "BeamBuilder.hxx"
#include "TestHelpers.hxx"
#include "TestMaterials.hxx"

#include "PrintSolution.hxx"

namespace
{
    //! Counts how many times each non-zero subdivision level fires.
    //! Index 0 corresponds to division level 1, index 1 to level 2, etc.
    class ObserveSimulationProgress : public Timesteps::TimestepObserver
    {
    public:
        void levelChanged(unsigned divisionLevel, unsigned) override
        {
            if(divisionLevel == 0)
            {
                return;
            }
            const size_t idx = divisionLevel - 1;
            if(idx >= m_SimulationCalls.size())
            {
                m_SimulationCalls.resize(idx + 1, 0u);
            }
            ++m_SimulationCalls[idx];
        }

        [[nodiscard]] const std::vector<unsigned> & calls() const
        {
            return m_SimulationCalls;
        }

    private:
        std::vector<unsigned> m_SimulationCalls;
    };
}   // namespace

class MultiDomain_BeamWall_Stucco : public testing::Test
{
protected:
    void TearDown() override
    {
        // Reset the SimulationProperties singleton so any flags this test
        // toggled (e.g. excludeWaterLiquidTransportation) do not bleed into
        // subsequent tests in the same binary.
        HygroThermFEM::SimulationProperties::Instance().reset();
    }
};

TEST_F(MultiDomain_BeamWall_Stucco, Stucco_99dot9_Percent)
{
    constexpr auto excludeWaterLiquidTransportation{false};
    constexpr auto excludeHeatOfEvaporation{true};
    constexpr auto excludeCapillaryConduction{true};
    constexpr auto excludeVaporDiffusionConduction{true};
    constexpr auto thermalConductivityMoistureAndTemperatureDependent{false};
    HygroThermFEM::SimulationProperties::Instance().setCalculationParameters(
      excludeWaterLiquidTransportation,
      excludeHeatOfEvaporation,
      excludeCapillaryConduction,
      excludeVaporDiffusionConduction,
      thermalConductivityMoistureAndTemperatureDependent);

    HygroThermFEM::MultiDomain multiDomain({.performThermal = false, .performMoisture = true});

    constexpr double initialTemperature = 30.0;
    constexpr double initialDomainHumidity = 0.999;
    constexpr double bcHumidity = 1.0;

    // Initial state: warm, near-saturation (the regime where the solver
    // currently struggles).
    constexpr HygroThermFEM::State initialState({.temperature = initialTemperature,
                                                 .humidity = initialDomainHumidity,
                                                 .pressure = 101325.0,
                                                 .liquidPercent = 1.0});

    // Register the stucco material so it is available to the elements.
    const auto & stucco = multiDomain.materials().createSolidMaterial(TestHelper::Stucco());

    ObserveSimulationProgress progressThermal;
    multiDomain.subscribeThermal(&progressThermal);

    ObserveSimulationProgress progressMoisture;
    multiDomain.subscribeMoisture(&progressMoisture);

    // Beam geometry:
    //   - 0.02 m of stucco, 3 elements wide
    //   - Total height 0.05 m, 1 element row
    TestHelper::BeamBuilder builder(multiDomain);
    builder.xStart(0.0)
      .height(0.05)
      .numElementsY(1)
      .state(initialState)
      .addSegment({.material = stucco.name(), .numElementsX = 3, .width = 0.02})
      .build();

    // Boundary conditions: identical to the initial domain state on both
    // edges. With no thermal or moisture gradient driving the problem the
    // solution should remain constant; any drift exposes solver issues.
    constexpr auto hc = 10.0;
    const HygroThermFEM::FixedBCHCCoefficients exteriorBc{
      /*airTemperature=*/initialTemperature, hc, /*airHumidity=*/bcHumidity};
    const HygroThermFEM::FixedBCHCCoefficients interiorBc{
      /*airTemperature=*/initialTemperature, hc, /*airHumidity=*/bcHumidity};

    for(auto [i1, i2] : builder.leftEdge())
    {
        multiDomain.createBC_FixedHc(i1, i2, exteriorBc);
    }
    for(auto [i1, i2] : builder.rightEdge())
    {
        multiDomain.createBC_FixedHc(i1, i2, interiorBc);
    }

    constexpr double dTime = 3600.0;
    constexpr int nSteps = 10;

    auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
    auto humidities = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);

    std::vector<std::vector<double>> temperatureSolution;
    std::vector<double> temperatureError;
    std::vector<std::vector<double>> waterContentSolution;
    std::vector<double> humidityError;


    for(int step = 0; step < nSteps; ++step)
    {
        const auto solution =
          multiDomain.transient(temperatures, humidities, dTime, static_cast<size_t>(step));

        temperatureSolution.push_back(solution.temperature);
        temperatureError.push_back(solution.temperatureError);
        waterContentSolution.push_back(solution.waterContent);
        humidityError.push_back(solution.humidityError);

        temperatures = solution.temperature;
        humidities = solution.humidity;
    }

    std::ofstream out("Stucco_99dot9.txt");
    TestHelper::printVector2D("waterContent", waterContentSolution, out);
    TestHelper::printVector("humidityError", humidityError, out);
    TestHelper::printVector("humidityIterations", progressMoisture.calls(), out);
}

TEST_F(MultiDomain_BeamWall_Stucco, Stucco_90_Percent)
{
    constexpr auto excludeWaterLiquidTransportation{false};
    constexpr auto excludeHeatOfEvaporation{true};
    constexpr auto excludeCapillaryConduction{true};
    constexpr auto excludeVaporDiffusionConduction{true};
    constexpr auto thermalConductivityMoistureAndTemperatureDependent{false};
    HygroThermFEM::SimulationProperties::Instance().setCalculationParameters(
      excludeWaterLiquidTransportation,
      excludeHeatOfEvaporation,
      excludeCapillaryConduction,
      excludeVaporDiffusionConduction,
      thermalConductivityMoistureAndTemperatureDependent);

    HygroThermFEM::MultiDomain multiDomain({.performThermal = false, .performMoisture = true});

    constexpr double initialTemperature = 30.0;
    constexpr double initialDomainHumidity = 0.9;
    constexpr double bcHumidity = 1.0;

    // Initial state: warm, near-saturation (the regime where the solver
    // currently struggles).
    constexpr HygroThermFEM::State initialState({.temperature = initialTemperature,
                                                 .humidity = initialDomainHumidity,
                                                 .pressure = 101325.0,
                                                 .liquidPercent = 1.0});

    // Register the stucco material so it is available to the elements.
    const auto & stucco = multiDomain.materials().createSolidMaterial(TestHelper::Stucco());

    ObserveSimulationProgress progressThermal;
    multiDomain.subscribeThermal(&progressThermal);

    ObserveSimulationProgress progressMoisture;
    multiDomain.subscribeMoisture(&progressMoisture);

    // Beam geometry:
    //   - 0.02 m of stucco, 3 elements wide
    //   - Total height 0.05 m, 1 element row
    TestHelper::BeamBuilder builder(multiDomain);
    builder.xStart(0.0)
      .height(0.05)
      .numElementsY(1)
      .state(initialState)
      .addSegment({.material = stucco.name(), .numElementsX = 3, .width = 0.02})
      .build();

    // Boundary conditions: identical to the initial domain state on both
    // edges. With no thermal or moisture gradient driving the problem the
    // solution should remain constant; any drift exposes solver issues.
    constexpr auto hc = 10.0;
    const HygroThermFEM::FixedBCHCCoefficients exteriorBc{
      /*airTemperature=*/initialTemperature, hc, /*airHumidity=*/bcHumidity};
    const HygroThermFEM::FixedBCHCCoefficients interiorBc{
      /*airTemperature=*/initialTemperature, hc, /*airHumidity=*/bcHumidity};

    for(auto [i1, i2] : builder.leftEdge())
    {
        multiDomain.createBC_FixedHc(i1, i2, exteriorBc);
    }
    for(auto [i1, i2] : builder.rightEdge())
    {
        multiDomain.createBC_FixedHc(i1, i2, interiorBc);
    }

    constexpr double dTime = 3600.0;
    constexpr int nSteps = 10;

    auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
    auto humidities = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);

    std::vector<std::vector<double>> temperatureSolution;
    std::vector<double> temperatureError;
    std::vector<std::vector<double>> waterContentSolution;
    std::vector<double> humidityError;


    for(int step = 0; step < nSteps; ++step)
    {
        const auto solution =
          multiDomain.transient(temperatures, humidities, dTime, static_cast<size_t>(step));

        temperatureSolution.push_back(solution.temperature);
        temperatureError.push_back(solution.temperatureError);
        waterContentSolution.push_back(solution.waterContent);
        humidityError.push_back(solution.humidityError);

        temperatures = solution.temperature;
        humidities = solution.humidity;
    }

    std::ofstream out("Stucco_90.txt");
    TestHelper::printVector2D("waterContent", waterContentSolution, out);
    TestHelper::printVector("humidityError", humidityError, out);
    TestHelper::printVector("humidityIterations", progressMoisture.calls(), out);
}

TEST_F(MultiDomain_BeamWall_Stucco, Stucco_50_Percent)
{
    constexpr auto excludeWaterLiquidTransportation{false};
    constexpr auto excludeHeatOfEvaporation{true};
    constexpr auto excludeCapillaryConduction{true};
    constexpr auto excludeVaporDiffusionConduction{true};
    constexpr auto thermalConductivityMoistureAndTemperatureDependent{false};
    HygroThermFEM::SimulationProperties::Instance().setCalculationParameters(
      excludeWaterLiquidTransportation,
      excludeHeatOfEvaporation,
      excludeCapillaryConduction,
      excludeVaporDiffusionConduction,
      thermalConductivityMoistureAndTemperatureDependent);

    HygroThermFEM::MultiDomain multiDomain({.performThermal = false, .performMoisture = true});

    constexpr double initialTemperature = 30.0;
    constexpr double initialDomainHumidity = 0.5;
    constexpr double bcHumidity = 1.0;

    // Initial state: warm, near-saturation (the regime where the solver
    // currently struggles).
    constexpr HygroThermFEM::State initialState({.temperature = initialTemperature,
                                                 .humidity = initialDomainHumidity,
                                                 .pressure = 101325.0,
                                                 .liquidPercent = 1.0});

    // Register the stucco material so it is available to the elements.
    const auto & stucco = multiDomain.materials().createSolidMaterial(TestHelper::Stucco());

    ObserveSimulationProgress progressThermal;
    multiDomain.subscribeThermal(&progressThermal);

    ObserveSimulationProgress progressMoisture;
    multiDomain.subscribeMoisture(&progressMoisture);

    // Beam geometry:
    //   - 0.02 m of stucco, 3 elements wide
    //   - Total height 0.05 m, 1 element row
    TestHelper::BeamBuilder builder(multiDomain);
    builder.xStart(0.0)
      .height(0.05)
      .numElementsY(1)
      .state(initialState)
      .addSegment({.material = stucco.name(), .numElementsX = 3, .width = 0.02})
      .build();

    // Boundary conditions: identical to the initial domain state on both
    // edges. With no thermal or moisture gradient driving the problem the
    // solution should remain constant; any drift exposes solver issues.
    constexpr auto hc = 10.0;
    const HygroThermFEM::FixedBCHCCoefficients exteriorBc{
      /*airTemperature=*/initialTemperature, hc, /*airHumidity=*/bcHumidity};
    const HygroThermFEM::FixedBCHCCoefficients interiorBc{
      /*airTemperature=*/initialTemperature, hc, /*airHumidity=*/bcHumidity};

    for(auto [i1, i2] : builder.leftEdge())
    {
        multiDomain.createBC_FixedHc(i1, i2, exteriorBc);
    }
    for(auto [i1, i2] : builder.rightEdge())
    {
        multiDomain.createBC_FixedHc(i1, i2, interiorBc);
    }

    constexpr double dTime = 3600.0;
    constexpr int nSteps = 10;

    auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
    auto humidities = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);

    std::vector<std::vector<double>> temperatureSolution;
    std::vector<double> temperatureError;
    std::vector<std::vector<double>> waterContentSolution;
    std::vector<double> humidityError;


    for(int step = 0; step < nSteps; ++step)
    {
        const auto solution =
          multiDomain.transient(temperatures, humidities, dTime, static_cast<size_t>(step));

        temperatureSolution.push_back(solution.temperature);
        temperatureError.push_back(solution.temperatureError);
        waterContentSolution.push_back(solution.waterContent);
        humidityError.push_back(solution.humidityError);

        temperatures = solution.temperature;
        humidities = solution.humidity;
    }

    std::ofstream out("Stucco_50.txt");
    TestHelper::printVector2D("waterContent", waterContentSolution, out);
    TestHelper::printVector("humidityError", humidityError, out);
    TestHelper::printVector("humidityIterations", progressMoisture.calls(), out);
}

TEST_F(MultiDomain_BeamWall_Stucco, Fiberglass_99dot99999_Percent)
{
    constexpr auto excludeWaterLiquidTransportation{false};
    constexpr auto excludeHeatOfEvaporation{true};
    constexpr auto excludeCapillaryConduction{true};
    constexpr auto excludeVaporDiffusionConduction{true};
    constexpr auto thermalConductivityMoistureAndTemperatureDependent{false};
    HygroThermFEM::SimulationProperties::Instance().setCalculationParameters(
      excludeWaterLiquidTransportation,
      excludeHeatOfEvaporation,
      excludeCapillaryConduction,
      excludeVaporDiffusionConduction,
      thermalConductivityMoistureAndTemperatureDependent);

    HygroThermFEM::MultiDomain multiDomain({.performThermal = false, .performMoisture = true});

    constexpr double initialTemperature = 30.0;
    constexpr double initialDomainHumidity = 0.9999999;
    constexpr double bcHumidity = 1.0;

    // Initial state: warm, near-saturation (the regime where the solver
    // currently struggles).
    constexpr HygroThermFEM::State initialState({.temperature = initialTemperature,
                                                 .humidity = initialDomainHumidity,
                                                 .pressure = 101325.0,
                                                 .liquidPercent = 1.0});

    // Register the fiberglass material so it is available to the elements.
    const auto & fiberglass =
      multiDomain.materials().createSolidMaterial(TestHelper::FiberglassBatts());

    ObserveSimulationProgress progressThermal;
    multiDomain.subscribeThermal(&progressThermal);

    ObserveSimulationProgress progressMoisture;
    multiDomain.subscribeMoisture(&progressMoisture);

    // Beam geometry:
    //   - 0.02 m of fiberglass, 3 elements wide
    //   - Total height 0.05 m, 1 element row
    TestHelper::BeamBuilder builder(multiDomain);
    builder.xStart(0.0)
      .height(0.05)
      .numElementsY(1)
      .state(initialState)
      .addSegment({.material = fiberglass.name(), .numElementsX = 3, .width = 0.02})
      .build();

    // Boundary conditions: identical to the initial domain state on both
    // edges. With no thermal or moisture gradient driving the problem the
    // solution should remain constant; any drift exposes solver issues.
    constexpr auto hc = 10.0;
    const HygroThermFEM::FixedBCHCCoefficients exteriorBc{
      /*airTemperature=*/initialTemperature, hc, /*airHumidity=*/bcHumidity};
    const HygroThermFEM::FixedBCHCCoefficients interiorBc{
      /*airTemperature=*/initialTemperature, hc, /*airHumidity=*/bcHumidity};

    for(auto [i1, i2] : builder.leftEdge())
    {
        multiDomain.createBC_FixedHc(i1, i2, exteriorBc);
    }
    for(auto [i1, i2] : builder.rightEdge())
    {
        multiDomain.createBC_FixedHc(i1, i2, interiorBc);
    }

    constexpr double dTime = 3600.0;
    constexpr int nSteps = 10;

    auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
    auto humidities = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);

    std::vector<std::vector<double>> temperatureSolution;
    std::vector<double> temperatureError;
    std::vector<std::vector<double>> waterContentSolution;
    std::vector<double> humidityError;


    for(int step = 0; step < nSteps; ++step)
    {
        const auto solution =
          multiDomain.transient(temperatures, humidities, dTime, static_cast<size_t>(step));

        temperatureSolution.push_back(solution.temperature);
        temperatureError.push_back(solution.temperatureError);
        waterContentSolution.push_back(solution.waterContent);
        humidityError.push_back(solution.humidityError);

        temperatures = solution.temperature;
        humidities = solution.humidity;
    }

    std::ofstream out("Fiberglass_99dot99999.txt");
    TestHelper::printVector2D("waterContent", waterContentSolution, out);
    TestHelper::printVector("humidityError", humidityError, out);
    TestHelper::printVector("humidityIterations", progressMoisture.calls(), out);
}
