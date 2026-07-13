#include <algorithm>
#include <cmath>
#include <numbers>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

#include "BeamBuilder.hxx"
#include "TestMaterials.hxx"

// D5 Tier 1: latent heat of fusion via the mass-conservative secant capacity
// (FusionSecantCapacity * condensed water), moisture immobilized. Mirrors the 1D
// reference solver's validation (hygrothermfem_python tests/test_freezing.py):
// the one-phase Stefan benchmark and the adiabatic two-phase equilibrium.
namespace
{
    //! Solves lam * exp(lam^2) * erf(lam) = Ste / sqrt(pi) by bisection.
    double stefanLambda(const double stefanNumber)
    {
        const double target = stefanNumber / std::sqrt(std::numbers::pi);
        double low = 1.0e-6;
        double high = 3.0;
        for(int iteration = 0; iteration < 200; ++iteration)
        {
            const double mid = 0.5 * (low + high);
            const double value = mid * std::exp(mid * mid) * std::erf(mid) - target;
            if(value > 0.0)
            {
                high = mid;
            }
            else
            {
                low = mid;
            }
        }
        return 0.5 * (low + high);
    }

    //! x where the bottom-row temperature profile crosses `level` (column-major nodes,
    //! numElementsY = 1 -> stride 2), linear interpolation.
    double frontPosition(const std::vector<double> & temperatures,
                         const double length,
                         const double level)
    {
        std::vector<double> profile;
        for(std::size_t idx = 0; idx < temperatures.size(); idx += 2)
        {
            profile.push_back(temperatures[idx]);
        }
        const double dx = length / static_cast<double>(profile.size() - 1);
        for(std::size_t idx = 1; idx < profile.size(); ++idx)
        {
            if(profile[idx - 1] < level && profile[idx] >= level)
            {
                const double fraction =
                  (level - profile[idx - 1]) / (profile[idx] - profile[idx - 1]);
                return (static_cast<double>(idx - 1) + fraction) * dx;
            }
        }
        return -1.0;
    }

    //! Immobilized-moisture setup: all four moisture-coupling terms excluded, fusion
    //! included. Must run before the domain geometry is built.
    void setFusionOnlyPhysics()
    {
        HygroThermFEM::SimulationProperties::Instance().setCalculationParameters(
          true, true, true, true, false);
        HygroThermFEM::SimulationProperties::Instance().setExcludeLatentHeatOfFusion(false);
    }

    void resetPhysics()
    {
        HygroThermFEM::SimulationProperties::Instance().resetCalculationParameters();
    }
}   // namespace

TEST(Freezing, StefanFrontPosition)
{
    // One-phase Stefan: a 0.5 m strip of Stucco at phi = 0.99, initially liquid AT the
    // fusion temperature (0 C); the left face drops to -10 C (Dirichlet via a huge film
    // coefficient). The frozen front must follow x = 2 lam sqrt(alpha t).
    setFusionOnlyPhysics();

    HygroThermFEM::MultiDomain multiDomain;
    multiDomain.performThermalSimulation(true);
    multiDomain.performMoistureSimulation(false);

    const auto & stucco = multiDomain.materials().createSolidMaterial(TestHelper::Stucco());

    constexpr double length{0.5};
    constexpr double initialHumidity{0.99};
    const HygroThermFEM::State initialState({.temperature = 0.0,
                                             .humidity = initialHumidity,
                                             .pressure = 101325.0,
                                             .liquidPercent = 1.0});

    TestHelper::BeamBuilder builder(multiDomain);
    builder.xStart(0.0)
      .height(0.05)
      .numElementsY(1)
      .state(initialState)
      .addSegment({.material = stucco.name(), .numElementsX = 200, .width = length})
      .build();

    constexpr double surfaceTemperature{-10.0};
    const HygroThermFEM::FixedBCHCCoefficients coldSurface{surfaceTemperature, 1.0e6, 0.5};
    builder.applyBC_FixedHc(TestHelper::BeamBuilder::Edge::Left, coldSurface);

    // Analytic front: frozen-region properties from the engine's capacitance convention
    // (rho_d + w)(c_d + w / rho_ice * cp_ice) with w = w(0.99) from the sorption curve.
    const double waterContent{95.0};   // Stucco sorption curve at phi = 0.99
    const double frozenCapacity =
      (1800.0 + waterContent) * (850.0 + waterContent / 916.7 * 2108.0);
    const double latentVolumetric = 333550.0 * waterContent;
    const double stefanNumber = frozenCapacity * 10.0 / latentVolumetric;
    const double lam = stefanLambda(stefanNumber);
    const double diffusivity = 0.85 / frozenCapacity;

    constexpr double dTime{600.0};
    constexpr unsigned numberOfSteps{144};

    auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
    auto humidities = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);
    std::vector<std::vector<double>> history{temperatures};
    for(unsigned step = 0; step < numberOfSteps; ++step)
    {
        const auto solution = multiDomain.transient(temperatures, humidities, dTime, step);
        temperatures = solution.temperature;
        humidities = solution.humidity;
        history.push_back(temperatures);
    }
    resetPhysics();

    for(const unsigned step : {48u, 96u, 144u})
    {
        const double analytic =
          2.0 * lam * std::sqrt(diffusivity * static_cast<double>(step) * dTime);
        const double numeric = frontPosition(history[step], length, -0.0005);
        EXPECT_NEAR(numeric, analytic, 0.05 * analytic) << "step " << step;
    }
}

TEST(Freezing, AdiabaticTwoPhaseEquilibrium)
{
    // Insulated strip, 60% of it liquid at +10 C and 40% frozen-cold at -10 C, high
    // water content. The sensible-only mean is clearly positive (~+2 C), but freezing
    // the excess enthalpy pins the equilibrium INSIDE the fusion ramp [-0.001, 0]:
    // reachable only if the fusion capacity actually absorbs the latent heat.
    setFusionOnlyPhysics();

    HygroThermFEM::MultiDomain multiDomain;
    multiDomain.performThermalSimulation(true);
    multiDomain.performMoistureSimulation(false);

    const auto & gypsum =
      multiDomain.materials().createSolidMaterial(TestHelper::GypsumBoardInterior());

    const HygroThermFEM::State coldState({.temperature = -10.0,
                                          .humidity = 0.999,
                                          .pressure = 101325.0,
                                          .liquidPercent = 0.0});

    TestHelper::BeamBuilder builder(multiDomain);
    builder.xStart(0.0)
      .height(0.05)
      .numElementsY(1)
      .state(coldState)
      .addSegment({.material = gypsum.name(), .numElementsX = 50, .width = 0.2})
      .build();

    // Rewrite the initial condition: first 40% of columns cold/frozen, rest warm/liquid.
    auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
    std::vector<double> liquidPercents(temperatures.size(), 0.0);
    for(std::size_t idx = 0; idx < temperatures.size(); ++idx)
    {
        const bool warmSide = (idx / 2) >= 20;   // column index (stride 2, numElementsY 1)
        temperatures[idx] = warmSide ? 10.0 : -10.0;
        liquidPercents[idx] = warmSide ? 1.0 : 0.0;
    }
    multiDomain.nodes().updateNodeTemperatures(temperatures, true);
    multiDomain.nodes().updateNodeLiquidPercents(liquidPercents, true);

    auto humidities = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);
    constexpr double dTime{300.0};
    constexpr unsigned numberOfSteps{2200};
    for(unsigned step = 0; step < numberOfSteps; ++step)
    {
        const auto solution = multiDomain.transient(temperatures, humidities, dTime, step);
        temperatures = solution.temperature;
        humidities = solution.humidity;
    }
    resetPhysics();

    const auto [minTemp, maxTemp] = std::minmax_element(temperatures.begin(), temperatures.end());
    EXPECT_GE(*minTemp, -0.0011);
    EXPECT_LE(*maxTemp, 0.0001);
}
