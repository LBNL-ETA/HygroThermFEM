#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

#include "BeamBuilder.hxx"
#include "TestMaterials.hxx"

// Steady state versus transient on one two-material beam, for each of THERM's moisture
// modelling options in turn.
//
// Beam: Stucco 76.3 mm (exterior) | Fiberglass Batts 80 mm (interior), 8 + 8 elements,
// one element tall, the wall of THERM's Stucco-Fiberglass pair. Constant NFRC films:
// -18 C / RH 0.5 / hc 26 on the stucco face, 21 C / RH 0.5 / hc 4.65 on the fiberglass face.
// Initial state 21 C / RH 0.5 everywhere. Same beam and BCs in every case.
//
// Five option sets: all off, then exactly one of water liquid transportation, heat of
// evaporation, capillary conduction, vapour diffusion (conduction) switched on.
//
// REFERENCE = the transient run, 8760 steps of 10 h (ten simulated years) under those
// constant films: by then the field has stopped changing and IS a steady solution of the
// same equations (after five years the stucco was still 0.005 RH short of rest). The transient results are hard-coded below, per x-column (the two nodes
// of a column carry the same value), so a reader sees the numbers. They were produced ONCE
// (2026-09-03, engine state at that date) and the transient tests that regenerate them are
// DISABLED_: rerun them only when the reference itself is in doubt
// (--gtest_also_run_disabled_tests, a few seconds per case).
//
// The steady-state tests run every time and chase those references.
namespace
{
    constexpr std::size_t columnCount{17u};

    struct Reference
    {
        std::vector<double> temperature;   //!< C, per x-column, exterior to interior
        std::vector<double> humidity;      //!< RH fraction, per x-column
    };

    struct ModellingCase
    {
        const char * name;
        bool liquidTransport;
        bool heatOfEvaporation;
        bool capillaryConduction;
        bool vaporDiffusionConduction;
        Reference transientReference;   //!< filled from the one-time transient run
    };

    void applyOptions(const ModellingCase & modellingCase)
    {
        HygroThermFEM::SimulationProperties::Instance().setCalculationParameters(
          !modellingCase.liquidTransport,
          !modellingCase.heatOfEvaporation,
          !modellingCase.capillaryConduction,
          !modellingCase.vaporDiffusionConduction,
          false);
    }

    void buildBeam(HygroThermFEM::MultiDomain & multiDomain)
    {
        const auto & stucco = multiDomain.materials().createSolidMaterial(TestHelper::Stucco());
        const auto & fiberglass = multiDomain.materials().createSolidMaterial(TestHelper::FiberglassBatts());

        TestHelper::BeamBuilder builder(multiDomain);
        builder.xStart(0.0)
          .height(0.0728)
          .numElementsY(1)
          .state({.temperature = 21.0, .humidity = 0.5, .pressure = 101325.0, .liquidPercent = 1.0})
          .addSegment({.material = stucco.name(), .numElementsX = 8, .width = 0.0763})
          .addSegment({.material = fiberglass.name(), .numElementsX = 8, .width = 0.08})
          .build();

        const HygroThermFEM::FixedBCHCCoefficients exterior{-18.0, 26.0, 0.5};
        const HygroThermFEM::FixedBCHCCoefficients interior{21.0, 4.65, 0.5};
        builder.applyBC_FixedHc(TestHelper::BeamBuilder::Edge::Left, exterior);
        builder.applyBC_FixedHc(TestHelper::BeamBuilder::Edge::Right, interior);
    }

    //! Every second node: the beam is one element tall, so nodes 2k and 2k+1 share a column.
    std::vector<double> perColumn(const std::vector<double> & nodal)
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
        Reference fields;
        bool converged;
        std::size_t passes;
    };

    SteadyRun runSteady(const ModellingCase & modellingCase)
    {
        applyOptions(modellingCase);
        HygroThermFEM::MultiDomain multiDomain;
        buildBeam(multiDomain);
        const auto solution = multiDomain.steadyState();
        HygroThermFEM::SimulationProperties::Instance().resetCalculationParameters();
        return {{perColumn(solution.temperature), perColumn(solution.humidity)},
                solution.converged,
                solution.steadyPasses};
    }

    //! The reference run: 8760 steps of 10 h (ten years) under the constant films.
    Reference runTransient(const ModellingCase & modellingCase)
    {
        constexpr std::size_t steps{8760u};
        constexpr double dTime{10.0 * 3600.0};

        applyOptions(modellingCase);
        HygroThermFEM::MultiDomain multiDomain;
        buildBeam(multiDomain);
        auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
        auto humidities = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);
        for(std::size_t step = 0u; step < steps; ++step)
        {
            const auto solution = multiDomain.transient(temperatures, humidities, dTime, step);
            temperatures = solution.temperature;
            humidities = solution.humidity;
        }
        HygroThermFEM::SimulationProperties::Instance().resetCalculationParameters();
        return {perColumn(temperatures), perColumn(humidities)};
    }

    double maxAbsDifference(const std::vector<double> & lhs, const std::vector<double> & rhs)
    {
        double result{0.0};
        for(std::size_t index = 0u; index < (std::min)(lhs.size(), rhs.size()); ++index)
        {
            result = (std::max)(result, std::abs(lhs[index] - rhs[index]));
        }
        return result;
    }

    void printSideBySide(const char * name, const Reference & reference, const Reference & candidate)
    {
        std::cout << "[Beam] " << name << ": column  T_ref  T_steady   RH_ref  RH_steady\n";
        for(std::size_t index = 0u; index < reference.humidity.size(); ++index)
        {
            std::cout << "[Beam] " << std::setw(6) << index + 1u << std::fixed << std::setprecision(2) << std::setw(8)
                      << reference.temperature[index] << std::setw(10) << candidate.temperature[index]
                      << std::setprecision(3) << std::setw(9) << reference.humidity[index] << std::setw(11)
                      << candidate.humidity[index] << "\n";
        }
        std::cout.unsetf(std::ios::fixed);
    }

    //! Prints a Reference as C++ initialisers, ready to paste into the case table.
    void printAsInitialisers(const char * name, const Reference & reference)
    {
        std::cout << "[Gen] " << name << "\n[Gen]   .temperature = {";
        for(const auto value : reference.temperature)
        {
            std::cout << std::setprecision(6) << value << ", ";
        }
        std::cout << "},\n[Gen]   .humidity = {";
        for(const auto value : reference.humidity)
        {
            std::cout << std::setprecision(6) << value << ", ";
        }
        std::cout << "},\n";
    }

    // ---------------------------------------------------------------------------------------
    // Reference tables: transient, 8760 x 10 h, generated 2026-09-03. Columns run from the
    // exterior stucco face (1) to the interior fiberglass face (17); column 9 is the
    // stucco/fiberglass interface.
    // ---------------------------------------------------------------------------------------
    const ModellingCase allOff{
      "all options off",
      false, false, false, false,
      {.temperature = {-17.3193, -17.1208, -16.9222, -16.7236, -16.5251, -16.3265, -16.1279, -15.9293, -15.7308,
                       -11.6151, -7.49953, -3.38392, 0.7317, 4.84732, 8.96293, 13.0785, 17.1942},
       .humidity = {0.473782, 0.547982, 0.619628, 0.688793, 0.755549, 0.819963, 0.882139, 0.942131, 1,
                    1, 1, 1, 0.949136, 0.866044, 0.771978, 0.678307, 0.590716}}};

    const ModellingCase liquidOnly{
      "water liquid transportation only",
      true, false, false, false,
      {.temperature = {-17.3193, -17.1208, -16.9222, -16.7236, -16.5251, -16.3265, -16.1279, -15.9293, -15.7308,
                       -11.6151, -7.49953, -3.38392, 0.7317, 4.84732, 8.96293, 13.0785, 17.1942},
       .humidity = {0.509261, 0.850784, 0.946495, 0.973749, 0.987747, 0.994606, 0.995885, 0.996491, 0.996951,
                    1, 1, 1, 0.949136, 0.866044, 0.771979, 0.678307, 0.590716}}};

    const ModellingCase evaporationOnly{
      "heat of evaporation only",
      false, true, false, false,
      {.temperature = {-17.1222, -16.8658, -16.6095, -16.3532, -16.0968, -15.8405, -15.5842, -15.3279, -15.0715,
                       -10.4289, -6.04372, -1.89431, 2.10169, 5.93982, 9.62075, 13.1485, 16.529},
       .humidity = {0.466389, 0.544422, 0.618971, 0.690162, 0.758117, 0.822953, 0.884812, 0.943793, 1,
                    1, 1, 0.973057, 0.908469, 0.831524, 0.754263, 0.682009, 0.616737}}};

    //! Capillary conduction is heat carried by the liquid flux; with liquid transport off
    //! there is no such flux, so the term is inactive and this case equals "all off". (Before
    //! 2026-09-03 the term ran on a fictitious flux built from the vapour-only humidity
    //! gradients and the transient blew up, -273 C mid-wall at 5 h steps.)
    const ModellingCase capillaryOnly{
      "capillary conduction only",
      false, false, true, false,
      {.temperature = {-17.3193, -17.1208, -16.9222, -16.7236, -16.5251, -16.3265, -16.1279, -15.9293, -15.7308,
                       -11.6151, -7.49953, -3.38392, 0.7317, 4.84732, 8.96293, 13.0785, 17.1942},
       .humidity = {0.473782, 0.547982, 0.619628, 0.688793, 0.755549, 0.819963, 0.882139, 0.942131, 1,
                    1, 1, 1, 0.949136, 0.866044, 0.771978, 0.678307, 0.590716}}};

    //! Capillary conduction is heat carried by the liquid flux, so it only exists together
    //! with liquid transport; this is the pair that exercises it.
    const ModellingCase liquidAndCapillary{
      "water liquid transportation + capillary conduction",
      true, false, true, false,
      {.temperature = {-17.3276, -17.1312, -16.9343, -16.737, -16.5394, -16.3416, -16.1434, -15.9449, -15.7463,
                       -11.6289, -7.51157, -3.39422, 0.723139, 4.84049, 8.95785, 13.0752, 17.1926},
       .humidity = {0.509575, 0.850806, 0.946488, 0.973736, 0.987731, 0.994605, 0.995884, 0.99649, 0.99695,
                    1, 1, 1, 0.949254, 0.866181, 0.772093, 0.678384, 0.590755}}};

    const ModellingCase vaporDiffusionOnly{
      "vapor diffusion only",
      false, false, false, true,
      {.temperature = {-17.3175, -17.1184, -16.9193, -16.7201, -16.521, -16.3219, -16.1228, -15.9237, -15.7246,
                       -11.5986, -7.4745, -3.35305, 0.765256, 4.8803, 8.99205, 13.1005, 17.2055},
       .humidity = {0.473711, 0.547947, 0.61962, 0.688804, 0.755571, 0.81999, 0.882163, 0.942146, 1,
                    1, 1, 1, 0.948276, 0.864917, 0.770949, 0.677573, 0.590364}}};

    //! Steady solve against the case's transient reference.
    void expectSteadyMatchesReference(const ModellingCase & modellingCase)
    {
        constexpr double temperatureTolerance{0.5};   // C
        constexpr double humidityTolerance{0.02};     // RH fraction (2 percentage points)

        const auto steady = runSteady(modellingCase);
        const auto & reference = modellingCase.transientReference;
        std::cout << "[Beam] " << modellingCase.name << ": steady converged = " << std::boolalpha << steady.converged
                  << ", passes = " << steady.passes << "\n";
        ASSERT_EQ(reference.humidity.size(), columnCount) << "reference table not generated for " << modellingCase.name;
        printSideBySide(modellingCase.name, reference, steady.fields);
        std::cout << "[Beam] max |dT| = " << maxAbsDifference(reference.temperature, steady.fields.temperature)
                  << ", max |dRH| = " << maxAbsDifference(reference.humidity, steady.fields.humidity) << "\n";

        EXPECT_TRUE(steady.converged) << modellingCase.name;
        EXPECT_LT(maxAbsDifference(reference.temperature, steady.fields.temperature), temperatureTolerance)
          << modellingCase.name;
        EXPECT_LT(maxAbsDifference(reference.humidity, steady.fields.humidity), humidityTolerance)
          << modellingCase.name;
    }

    //! Regenerates the transient reference and compares it with the hard-coded table.
    void expectTransientMatchesTable(const ModellingCase & modellingCase)
    {
        constexpr double tolerance{1e-4};
        const auto transient = runTransient(modellingCase);
        printAsInitialisers(modellingCase.name, transient);
        ASSERT_EQ(modellingCase.transientReference.humidity.size(), columnCount)
          << "reference table not generated for " << modellingCase.name;
        EXPECT_LT(maxAbsDifference(modellingCase.transientReference.temperature, transient.temperature), tolerance);
        EXPECT_LT(maxAbsDifference(modellingCase.transientReference.humidity, transient.humidity), tolerance);
    }
}   // namespace

// Steady-state cases: run every time, chased against the transient references.
TEST(SteadyVsTransient_StuccoFiberglassBeam, Steady_AllOff)
{
    expectSteadyMatchesReference(allOff);
}

TEST(SteadyVsTransient_StuccoFiberglassBeam, Steady_LiquidTransportOnly)
{
    expectSteadyMatchesReference(liquidOnly);
}

TEST(SteadyVsTransient_StuccoFiberglassBeam, Steady_HeatOfEvaporationOnly)
{
    expectSteadyMatchesReference(evaporationOnly);
}

TEST(SteadyVsTransient_StuccoFiberglassBeam, Steady_CapillaryConductionOnly)
{
    expectSteadyMatchesReference(capillaryOnly);
}

TEST(SteadyVsTransient_StuccoFiberglassBeam, Steady_LiquidTransportAndCapillaryConduction)
{
    expectSteadyMatchesReference(liquidAndCapillary);
}

TEST(SteadyVsTransient_StuccoFiberglassBeam, Steady_VaporDiffusionOnly)
{
    expectSteadyMatchesReference(vaporDiffusionOnly);
}

// Transient references: run once; rerun only when the reference is in doubt.
TEST(SteadyVsTransient_StuccoFiberglassBeam, DISABLED_Transient_AllOff)
{
    expectTransientMatchesTable(allOff);
}

TEST(SteadyVsTransient_StuccoFiberglassBeam, DISABLED_Transient_LiquidTransportOnly)
{
    expectTransientMatchesTable(liquidOnly);
}

TEST(SteadyVsTransient_StuccoFiberglassBeam, DISABLED_Transient_HeatOfEvaporationOnly)
{
    expectTransientMatchesTable(evaporationOnly);
}

TEST(SteadyVsTransient_StuccoFiberglassBeam, DISABLED_Transient_CapillaryConductionOnly)
{
    expectTransientMatchesTable(capillaryOnly);
}

TEST(SteadyVsTransient_StuccoFiberglassBeam, DISABLED_Transient_VaporDiffusionOnly)
{
    expectTransientMatchesTable(vaporDiffusionOnly);
}

TEST(SteadyVsTransient_StuccoFiberglassBeam, DISABLED_Transient_LiquidTransportAndCapillaryConduction)
{
    expectTransientMatchesTable(liquidAndCapillary);
}
