#include <algorithm>
#include <cmath>
#include <vector>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

#include "BeamBuilder.hxx"
#include "TestMaterials.hxx"

// D5 Tier 1 conservation guard: repeated freeze--thaw cycling through the fusion ramp.
// The Freezing tests only ever freeze monotonically; here every near-surface node
// crosses the ramp in BOTH directions several times, which is where an asymmetric or
// non-telescoping latent booking would hide. Mirrors the 1D reference solver's
// freeze--thaw enthalpy-conservation test (hygrothermfem_python tests/test_freezing.py).
//
// Checks, in order: freezing actually happens each cold phase; the slab returns to its
// initial uniform state after the final warm hold (no hysteresis -- thaw repays exactly
// what freezing released); and the film-flux integral matches the enthalpy inventory
// (the energy audit that exposed the false-convergence leak, made permanent).
namespace
{
    constexpr double slabLength{0.1};
    constexpr double slabHeight{0.05};
    constexpr std::size_t elementsX{50};
    constexpr double dTime{600.0};
    constexpr double filmCoefficient{1000.0};
    constexpr double initialTemperature{5.0};

    constexpr unsigned coldSteps{36};
    constexpr unsigned warmSteps{36};
    constexpr unsigned numberOfCycles{3};
    constexpr unsigned holdSteps{144};
    constexpr unsigned totalSteps{numberOfCycles * (coldSteps + warmSteps) + holdSteps};

    //! Air temperature series: numberOfCycles of (-10 C, +10 C), then a +5 C hold long
    //! enough for the slab to relax back to its uniform initial state.
    std::vector<double> airTemperatureSeries()
    {
        std::vector<double> series;
        for(unsigned cycle = 0; cycle < numberOfCycles; ++cycle)
        {
            series.insert(series.end(), coldSteps, -10.0);
            series.insert(series.end(), warmSteps, 10.0);
        }
        series.insert(series.end(), holdSteps, initialTemperature);
        return series;
    }

    std::vector<HygroThermFEM::FixedBCHCCoefficients> filmSeries()
    {
        std::vector<HygroThermFEM::FixedBCHCCoefficients> series;
        for(const double airTemperature : airTemperatureSeries())
        {
            series.push_back({airTemperature, filmCoefficient, 0.8});
        }
        return series;
    }

    //! Nodal volumes of the regular 2-row strip mesh (per metre of depth): each column
    //! owns dx * height (half at the ends), split evenly between its two row nodes.
    std::vector<double> nodalVolumes()
    {
        const double dx = slabLength / static_cast<double>(elementsX);
        std::vector<double> volumes;
        for(std::size_t column = 0; column <= elementsX; ++column)
        {
            const bool endColumn = (column == 0) || (column == elementsX);
            const double columnWidth = endColumn ? dx / 2.0 : dx;
            volumes.push_back(columnWidth * slabHeight / 2.0);
            volumes.push_back(columnWidth * slabHeight / 2.0);
        }
        return volumes;
    }

    //! Nodal enthalpy per unit volume, zero for liquid at 0 C. Mirrors the engine's
    //! capacitance convention for Stucco at phi = 0.99 (w = 95 kg/m3): the sensible
    //! (rho_d + w)(c_d + w/rho * cp) split between ice and water, plus the fusion jump
    //! L_f * w across the ramp [-0.1, 0] C.
    double nodalEnthalpy(const double temperature)
    {
        constexpr double waterContent{95.0};
        constexpr double frozenCapacity =
          (1800.0 + waterContent) * (850.0 + waterContent / 916.7 * 2108.0);
        constexpr double liquidCapacity =
          (1800.0 + waterContent) * (850.0 + waterContent / 1000.0 * 4187.0);
        constexpr double latentVolumetric = 333550.0 * waterContent;
        constexpr double rampBottom{-0.1};

        if(temperature >= 0.0)
        {
            return liquidCapacity * temperature;
        }
        if(temperature <= rampBottom)
        {
            return frozenCapacity * (temperature - rampBottom)
                   + 0.5 * (frozenCapacity + liquidCapacity) * rampBottom - latentVolumetric;
        }
        const double liquidFraction = (temperature - rampBottom) / (-rampBottom);
        return 0.5 * (frozenCapacity + liquidCapacity) * temperature
               + latentVolumetric * (liquidFraction - 1.0);
    }

    double slabEnthalpy(const std::vector<double> & temperatures,
                        const std::vector<double> & volumes)
    {
        double total{0.0};
        for(std::size_t idx = 0; idx < temperatures.size(); ++idx)
        {
            total += volumes[idx] * nodalEnthalpy(temperatures[idx]);
        }
        return total;
    }

    //! Film flux booked by a backward-Euler step: h * (T_air - T_surface_end) over the
    //! left edge (two rows, each owning half the edge height).
    double stepFilmEnergy(const std::vector<double> & temperatures, const double airTemperature)
    {
        const double edgeWeight = slabHeight / 2.0;
        const double fluxRow0 = filmCoefficient * (airTemperature - temperatures[0]);
        const double fluxRow1 = filmCoefficient * (airTemperature - temperatures[1]);
        return (fluxRow0 + fluxRow1) * edgeWeight * dTime;
    }
}   // namespace

TEST(FreezeThawCycling, EnthalpyConservationAndNoHysteresis)
{
    HygroThermFEM::SimulationProperties::Instance().setCalculationParameters(
      true, true, true, true, false);
    HygroThermFEM::SimulationProperties::Instance().setExcludeLatentHeatOfFusion(false);

    HygroThermFEM::MultiDomain multiDomain;
    multiDomain.performThermalSimulation(true);
    multiDomain.performMoistureSimulation(false);

    const auto & stucco = multiDomain.materials().createSolidMaterial(TestHelper::Stucco());
    const HygroThermFEM::State initialState({.temperature = initialTemperature,
                                             .humidity = 0.99,
                                             .pressure = 101325.0,
                                             .liquidPercent = 1.0});

    TestHelper::BeamBuilder builder(multiDomain);
    builder.xStart(0.0)
      .height(slabHeight)
      .numElementsY(1)
      .state(initialState)
      .addSegment({.material = stucco.name(), .numElementsX = elementsX, .width = slabLength})
      .build();
    builder.applyBC_FixedHc(TestHelper::BeamBuilder::Edge::Left, filmSeries());

    const auto airSeries = airTemperatureSeries();
    const auto volumes = nodalVolumes();

    auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
    auto humidities = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);
    ASSERT_EQ(volumes.size(), temperatures.size());

    const double initialEnthalpy = slabEnthalpy(temperatures, volumes);
    double netFilmEnergy{0.0};
    double grossFilmEnergy{0.0};
    for(unsigned step = 0; step < totalSteps; ++step)
    {
        const auto solution = multiDomain.transient(temperatures, humidities, dTime, step);
        temperatures = solution.temperature;
        humidities = solution.humidity;

        const double stepEnergy = stepFilmEnergy(temperatures, airSeries[step]);
        netFilmEnergy += stepEnergy;
        grossFilmEnergy += std::abs(stepEnergy);

        const bool insideCycles = step < numberOfCycles * (coldSteps + warmSteps);
        if(insideCycles && (step + 1) % (coldSteps + warmSteps) == coldSteps)
        {
            // End of a cold phase: the surface region must actually have frozen.
            const auto minTemp = *std::min_element(temperatures.begin(), temperatures.end());
            EXPECT_LT(minTemp, -2.0) << "cold phase ending at step " << step + 1;
        }
    }
    HygroThermFEM::SimulationProperties::Instance().resetCalculationParameters();

    // No hysteresis: after the long hold at the initial air temperature every node is
    // back at the initial state -- thaw absorbed exactly the latent heat freezing
    // released, nothing stayed frozen and no spurious energy accumulated.
    const auto [minTemp, maxTemp] = std::minmax_element(temperatures.begin(), temperatures.end());
    EXPECT_NEAR(*minTemp, initialTemperature, 0.05);
    EXPECT_NEAR(*maxTemp, initialTemperature, 0.05);

    // Energy closure: the film-flux integral must match the enthalpy inventory change.
    // This is the audit that exposed the false-convergence leak (which lost ~6 % of the
    // gross throughput); the bound also absorbs the step-end flux approximation when
    // the solver substeps internally.
    const double finalEnthalpy = slabEnthalpy(temperatures, volumes);
    const double imbalance = (finalEnthalpy - initialEnthalpy) - netFilmEnergy;
    EXPECT_LT(std::abs(imbalance), 0.02 * grossFilmEnergy);
}
