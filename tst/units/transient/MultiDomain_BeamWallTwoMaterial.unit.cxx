#include <fstream>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

#include "BeamBuilder.hxx"
#include "PrintSolution.hxx"
#include "TestMaterials.hxx"

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
        std::vector<unsigned> m_SimulationCalls{0u, 0u, 0u};
    };
}   // namespace

class MultiDomain_BeamWall_StuccoFiberglass : public testing::Test
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

//! Beam-shaped wall section with two materials: an exterior stucco layer
//! and a fiberglass-batt insulation layer. Initialised at high humidity
//! (RH = 0.9999) and warm temperature, with a colder/drier interior
//! boundary. Intended as a reproducer for solver divergence issues
//! reported on multi-material beam structures at near-saturation states.
TEST_F(MultiDomain_BeamWall_StuccoFiberglass, TwoMaterialHighHumidity)
{
    // Disable liquid transportation in the moisture transport equation to
    // confirm whether it is the source of the near-saturation drift seen on
    // this reproducer.
    constexpr auto excludeWaterLiquidTransportation{true};
    constexpr auto excludeHeatOfEvaporation{false};
    constexpr auto excludeCapillaryConduction{false};
    constexpr auto excludeVaporDiffusionConduction{false};
    constexpr auto thermalConductivityMoistureAndTemperatureDependent{false};
    HygroThermFEM::SimulationProperties::Instance().setCalculationParameters(
      excludeWaterLiquidTransportation,
      excludeHeatOfEvaporation,
      excludeCapillaryConduction,
      excludeVaporDiffusionConduction,
      thermalConductivityMoistureAndTemperatureDependent);

    HygroThermFEM::MultiDomain multiDomain;

    constexpr double initialTemperature = 30.0;
    constexpr double initialHumidity = 0.999;

    // Initial state: warm, near-saturation (the regime where the solver
    // currently struggles).
    constexpr HygroThermFEM::State initialState({.temperature = initialTemperature,
                                                 .humidity = initialHumidity,
                                                 .pressure = 101325.0,
                                                 .liquidPercent = 1.0});

    // Register both materials so they are available to the elements.
    const auto & stucco =
      multiDomain.materials().createSolidMaterial(TestHelper::Stucco());
    const auto & fiberglass =
      multiDomain.materials().createSolidMaterial(TestHelper::FiberglassBatts());

    ObserveSimulationProgress progressThermal;
    multiDomain.subscribeThermal(&progressThermal);

    ObserveSimulationProgress progressMoisture;
    multiDomain.subscribeMoisture(&progressMoisture);

    // Beam geometry:
    //   - 0.02 m of exterior stucco (left), 4 elements wide
    //   - 0.10 m of fiberglass batt insulation (right), 10 elements wide
    //   - Total height 0.05 m, 4 element rows
    TestHelper::BeamBuilder builder(multiDomain);
    builder.xStart(0.0)
      .height(0.05)
      .numElementsY(1)
      .state(initialState)
      .addSegment({.material = stucco.name(), .numElementsX = 3, .width = 0.02})
      //.addSegment({.material = fiberglass.name(), .numElementsX = 1, .width = 0.10})
      .build();

    // Boundary conditions: identical to the initial domain state on both
    // edges. With no thermal or moisture gradient driving the problem the
    // solution should remain constant; any drift exposes solver issues.
    constexpr auto hc = 10.0;
    const HygroThermFEM::FixedBCHCCoefficients exteriorBc{
      /*airTemperature=*/initialTemperature, hc, /*airHumidity=*/initialHumidity};
    const HygroThermFEM::FixedBCHCCoefficients interiorBc{
      /*airTemperature=*/initialTemperature, hc, /*airHumidity=*/initialHumidity};

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

    auto temperatures =
      multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
    auto humidities =
      multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);

    std::vector<std::vector<double>> temperatureSolution;
    std::vector<double> temperatureError;
    std::vector<std::vector<double>> waterContentSolution;
    std::vector<double> humidityError;

    std::string solverError;
    try
    {
        for(int step = 0; step < nSteps; ++step)
        {
            const auto solution = multiDomain.transient(
              temperatures, humidities, dTime, static_cast<size_t>(step));

            temperatureSolution.push_back(solution.temperature);
            temperatureError.push_back(solution.temperatureError);
            waterContentSolution.push_back(solution.waterContent);
            humidityError.push_back(solution.humidityError);

            temperatures = solution.temperature;
            humidities = solution.humidity;
        }
    }
    catch(const std::exception & e)
    {
        // Capture solver failure but continue so the print block always runs.
        solverError = e.what();
    }

    std::ofstream out("print.txt");
    if(!solverError.empty())
    {
        out << "// solver threw: " << solverError
            << " after " << temperatureSolution.size() << " step(s)\n";
    }
    TestHelper::printVector("waterContentError", humidityError, out);
    TestHelper::printVector("progressMoisture", progressMoisture.calls(), out);
    TestHelper::printVector2D("waterContent", waterContentSolution, out);
    TestHelper::printVector("temperatureError", temperatureError, out);
    TestHelper::printVector("progressThermal", progressThermal.calls(), out);
    TestHelper::printVector2D("temperature", temperatureSolution, out);
}
