#include <algorithm>
#include <cmath>
#include <numbers>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

#include "BeamBuilder.hxx"
#include "TestMaterials.hxx"

// D5 Tier 1: latent heat of fusion via the mass-conservative secant capacity
// (FusionSecantCapacity * condensed water), moisture immobilized. Validated the way
// the enthalpy method usually is: the one-phase Stefan benchmark and the adiabatic
// two-phase equilibrium.
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

// Anchor the helper above to EXTERNALLY published values before trusting it as the
// reference for everything below. stefanLambda() is an independent implementation
// of that root, so it needs its own anchor: another implementation being correct
// says nothing about this one. A dropped sqrt(pi) or a bracket that excludes the
// root would make StefanFrontPosition validate the engine against a wrong front
// and still pass.
//
// The published values are cited in full at the assertions below.
TEST(Freezing, StefanRootMatchesPublishedValues)
{
    // Gobin & Le Quere, Computer Assisted Mechanics and Engineering Sciences 7(3),
    // 289-306 (2000), p.305. Published expressly as the yardstick for judging code
    // accuracy in the pure-conduction limit of a 13-team comparison exercise.
    EXPECT_NEAR(stefanLambda(0.01), 0.0705932, 1.0e-7);
    EXPECT_NEAR(stefanLambda(0.1), 0.2200163, 1.0e-7);

    // Lunardini, Heat Conduction with Freezing or Thawing, CRREL Monograph 88-1
    // (1988), Table 2.2 p.32, zero-superheat column. Printed values are truncated
    // rather than rounded, so 1e-3 is the meaningful tolerance.
    EXPECT_NEAR(stefanLambda(1.0), 0.6201, 1.0e-3);

    // Alexiades & Solomon, Mathematical Modeling of Melting and Freezing Processes
    // (1993), section 2.1.E worked example: Ste = 0.314 -> lam = 0.3777 (4 d.p.).
    EXPECT_NEAR(stefanLambda(0.314), 0.3777, 1.0e-4);
}

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
        const double numeric = frontPosition(history[step], length, -0.05);
        EXPECT_NEAR(numeric, analytic, 0.05 * analytic) << "step " << step;
    }
}

TEST(Freezing, IceContentReportedBehindFront)
{
    // The Stefan tests validate the thermal front, which depends only on the TOTAL
    // condensed water (liquid + ice). This one pins the reported split: behind the
    // front the condensed water must be reported as ice, ahead of it as liquid. A
    // stalled liquid fraction leaves the front intact and the ice field identically
    // zero, which is exactly what the tests above cannot see.
    setFusionOnlyPhysics();

    HygroThermFEM::MultiDomain multiDomain;
    multiDomain.performThermalSimulation(true);
    multiDomain.performMoistureSimulation(false);

    const auto & stucco = multiDomain.materials().createSolidMaterial(TestHelper::Stucco());

    constexpr double length{0.5};
    const HygroThermFEM::State initialState({.temperature = 0.0,
                                             .humidity = 0.99,
                                             .pressure = 101325.0,
                                             .liquidPercent = 1.0});

    TestHelper::BeamBuilder builder(multiDomain);
    builder.xStart(0.0)
      .height(0.05)
      .numElementsY(1)
      .state(initialState)
      .addSegment({.material = stucco.name(), .numElementsX = 200, .width = length})
      .build();

    const HygroThermFEM::FixedBCHCCoefficients coldSurface{-10.0, 1.0e6, 0.5};
    builder.applyBC_FixedHc(TestHelper::BeamBuilder::Edge::Left, coldSurface);

    constexpr double dTime{600.0};
    constexpr unsigned numberOfSteps{48};

    auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
    auto humidities = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);
    HygroThermFEM::Solution solution = multiDomain.currentStateSolution();
    for(unsigned step = 0; step < numberOfSteps; ++step)
    {
        solution = multiDomain.transient(temperatures, humidities, dTime, step);
        temperatures = solution.temperature;
        humidities = solution.humidity;
    }
    resetPhysics();

    ASSERT_EQ(solution.iceContent.size(), solution.temperature.size());

    // Column-major nodes, one element row: node 0 sits on the cold face, the last node
    // on the far (still liquid) end.
    const std::size_t coldFace{0u};
    const std::size_t farEnd{solution.temperature.size() - 1u};
    EXPECT_LT(solution.temperature[coldFace], -5.0);
    EXPECT_GT(solution.temperature[farEnd], -0.1);

    const double condensedAtColdFace =
      solution.waterContent[coldFace] - solution.vaporContent[coldFace];
    EXPECT_GT(condensedAtColdFace, 50.0);
    EXPECT_NEAR(solution.iceContent[coldFace], condensedAtColdFace, 1e-6);
    EXPECT_NEAR(solution.liquidWaterContent[coldFace], 0.0, 1e-6);

    const double condensedAtFarEnd = solution.waterContent[farEnd] - solution.vaporContent[farEnd];
    EXPECT_NEAR(solution.iceContent[farEnd], 0.0, 1e-6);
    EXPECT_NEAR(solution.liquidWaterContent[farEnd], condensedAtFarEnd, 1e-6);
}

TEST(Freezing, IceContentReportedWithMoistureCoupling)
{
    // THERM's transient configuration: both domains solved, every moisture coupling term
    // on, moisture/temperature-dependent conductivity on, fusion on. A saturated stucco
    // strip between a -10 C face and a warm face must report ice at the cold face.
    HygroThermFEM::SimulationProperties::Instance().setCalculationParameters(
      false, false, false, false, true);
    HygroThermFEM::SimulationProperties::Instance().setExcludeLatentHeatOfFusion(false);

    HygroThermFEM::MultiDomain multiDomain({.performThermal = true, .performMoisture = true});

    const auto & stucco = multiDomain.materials().createSolidMaterial(TestHelper::Stucco());

    const HygroThermFEM::State initialState({.temperature = 21.0,
                                             .humidity = 0.99,
                                             .pressure = 101325.0,
                                             .liquidPercent = 1.0});

    TestHelper::BeamBuilder builder(multiDomain);
    builder.xStart(0.0)
      .height(0.05)
      .numElementsY(1)
      .state(initialState)
      .addSegment({.material = stucco.name(), .numElementsX = 40, .width = 0.2})
      .build();

    // Strong cold film on a 0.2 m strip: the steady cold-face temperature is about
    // -16 C, so the face is well inside the frozen range long before 24 h.
    const HygroThermFEM::FixedBCHCCoefficients coldSurface{-20.0, 25.0, 0.1};
    const HygroThermFEM::FixedBCHCCoefficients warmSurface{20.0, 8.0, 0.3};
    builder.applyBC_FixedHc(TestHelper::BeamBuilder::Edge::Left, coldSurface);
    builder.applyBC_FixedHc(TestHelper::BeamBuilder::Edge::Right, warmSurface);

    constexpr double dTime{3600.0};
    constexpr unsigned numberOfSteps{24};

    auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
    auto humidities = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);
    HygroThermFEM::Solution solution = multiDomain.currentStateSolution();
    for(unsigned step = 0; step < numberOfSteps; ++step)
    {
        solution = multiDomain.transient(temperatures, humidities, dTime, step);
        temperatures = solution.temperature;
        humidities = solution.humidity;
    }
    resetPhysics();

    const std::size_t coldFace{0u};
    ASSERT_LT(solution.temperature[coldFace], -0.1) << "cold face never froze";

    // The face dries toward the 10 % ambient over the day, so only a few kg/m3 remain
    // condensed there; the point is that whatever remains is reported frozen.
    const double condensedAtColdFace =
      solution.waterContent[coldFace] - solution.vaporContent[coldFace];
    EXPECT_GT(condensedAtColdFace, 1.0);
    EXPECT_NEAR(solution.iceContent[coldFace], condensedAtColdFace, 1e-6);
    EXPECT_NEAR(solution.liquidWaterContent[coldFace], 0.0, 1e-6);
}

TEST(Freezing, AdiabaticTwoPhaseEquilibrium)
{
    // Insulated strip, 60% of it liquid at +10 C and 40% frozen-cold at -10 C, high
    // water content. The sensible-only mean is clearly positive (~+2 C), but freezing
    // the excess enthalpy pins the equilibrium INSIDE the fusion ramp [-0.1, 0]:
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
    EXPECT_GE(*minTemp, -0.101);
    EXPECT_LE(*maxTemp, 0.0001);
}
