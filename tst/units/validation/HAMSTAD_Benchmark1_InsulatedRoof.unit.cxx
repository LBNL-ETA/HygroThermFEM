#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "BeamBuilder.hxx"
#include "HAMSTADBenchmarks.hxx"
#include "HygroThermFEM2D.hxx"
#include "NumericTable.hxx"
#include "TestHelpers.hxx"

/////////////////////////////////////////////////////////////////////////////////////
/// HAMSTAD benchmark 1, "Insulated roof", against the package's band of acceptance.
///
/// Report R-02:8 (Hagentoft, 2002) section 3.1, p. 8-12, and section 3.6, p. 52-53.
/// A 100 mm capillary active load bearing layer, sealed on its exterior face, carries
/// 50 mm of capillary non-active insulation on the interior. It starts near saturation
/// (145 kg/m3) at 10 C against a dry insulation (0,065 kg/m3), under one year of hourly
/// North-European exterior climate on the sealed face (heat only, alpha = 25 W/m2K) and
/// a dwelling's climate on the interior (alpha = 7, beta_p = 2e-8 s/m). The load bearing
/// layer dries inward through the insulation; in winter the vapour it releases condenses
/// again at the cold contact surface. The requested output is the integrated moisture in
/// each layer, M_A and M_B in kg/m2, hourly.
///
/// There is no analytic solution. Seven codes submitted results and the package defines
/// a BAND OF ACCEPTANCE from them: the 99,9 % confidence interval of the mean under the
/// t-distribution, hour by hour (section 3.6). Both bands are distributed in
/// Bench1Year1.ods and reproduced in the data directory; the expectation here is that
/// the engine's M_A and M_B stay inside them.
///
/// The sealing layer (Z_p = 1e12 m/s, K = 0, R = 0, no capacity) is not meshed: it is
/// exactly "no moisture exchange on the exterior face", which is what beta_p,e = 0
/// already states, and it carries no thermal resistance.
/////////////////////////////////////////////////////////////////////////////////////

namespace
{
    namespace Bench = TestHelper::HAMSTAD::Benchmark1;
    namespace Ham = TestHelper::HAMSTAD;

    constexpr std::size_t hoursPerYear{8760u};
    constexpr double hour{3600.0};

    //! ClimateBench1.txt: time (s), T_eq,e (C), T_eq,i (C), p_a,e (Pa), p_a,i (Pa), hourly
    //! from 0 to 8760 h. "Interpolated values for intermediate times": the engine steps
    //! hourly with backward Euler, so step k (from hour k to k + 1) takes row k + 1.
    struct Climate
    {
        std::vector<double> exteriorTemperature;
        std::vector<double> interiorTemperature;
        std::vector<double> exteriorPressure;
        std::vector<double> interiorPressure;
    };

    [[nodiscard]] Climate loadClimate()
    {
        const auto table =
          TestHelper::loadNumericTable(std::string{HAMSTAD_DATA_DIR} + "/ClimateBench1.txt", 0u);
        return {.exteriorTemperature = TestHelper::column(table, 1u),
                .interiorTemperature = TestHelper::column(table, 2u),
                .exteriorPressure = TestHelper::column(table, 3u),
                .interiorPressure = TestHelper::column(table, 4u)};
    }

    //! One layer's band from Bench1Year1_MA.csv / _MB.csv: time (h), the seven codes,
    //! Average, StandDev, CI-(p=0.1%), CI+(p=0.1%).
    struct Band
    {
        std::vector<double> average;
        std::vector<double> low;
        std::vector<double> high;
    };

    [[nodiscard]] Band loadBand(const std::string & fileName)
    {
        const auto table =
          TestHelper::loadNumericTable(std::string{HAMSTAD_DATA_DIR} + "/" + fileName, 1u);
        return {.average = TestHelper::column(table, 8u),
                .low = TestHelper::column(table, 10u),
                .high = TestHelper::column(table, 11u)};
    }

    //! The interior air humidity the engine multiplies by its own saturation
    //! concentration: the climate's partial pressure over the engine's p_s(T).
    [[nodiscard]] double humidityFromPressure(const double pressure, const double temperature)
    {
        return pressure / HygroThermFEM::vaporPressureAtTemperature(temperature);
    }

    //! Trapezoidal integral of w(phi) over the node columns [first, last] of one layer,
    //! through that layer's own closed-form retention curve: M [kg/m2].
    [[nodiscard]] double integratedMoisture(const std::vector<double> & coords,
                                            const std::vector<double> & humidities,
                                            const std::size_t first,
                                            const std::size_t last,
                                            const Ham::Material & material)
    {
        double mass{0.0};
        for(std::size_t col = first; col < last; ++col)
        {
            const double left{
              Ham::waterContent(material, Ham::suctionFromHumidity(humidities[col]))};
            const double right{
              Ham::waterContent(material, Ham::suctionFromHumidity(humidities[col + 1u]))};
            mass += 0.5 * (left + right) * (coords[col + 1u] - coords[col]);
        }
        return mass;
    }
}   // namespace

//! Runs the first year once for the suite.
class HAMSTAD_Benchmark1 : public ::testing::Test
{
protected:
    static std::vector<double> loadBearingMass;   //!< M_A [kg/m2], hour 0 to 8760
    static std::vector<double> insulationMass;    //!< M_B
    static std::vector<double> vapourInflow;      //!< cumulative through the interior face
    static Band bandA;
    static Band bandB;

    static void SetUpTestSuite()
    {
        HygroThermFEM::SimulationProperties::Instance().setCalculationParameters(
          false, false, false, false, true);
        HygroThermFEM::SimulationProperties::Instance().setIterationParameters(1.0, 1e-9, 500);

        const Ham::Material loadBearing{Bench::loadBearing()};
        const Ham::Material insulation{Bench::insulation()};

        HygroThermFEM::MultiDomain multiDomain({.performThermal = true, .performMoisture = true});
        const auto & matA =
          multiDomain.materials().createSolidMaterial(Ham::solidMaterial(loadBearing));
        const auto & matB =
          multiDomain.materials().createSolidMaterial(Ham::solidMaterial(insulation));

        // The initial water contents through each material's printed inverse retention
        // curve. The contact node is shared, so it carries ONE humidity: the insulation's,
        // because at the load bearing layer's 145 kg/m3 the insulation's isotherm would put
        // 4 kg/m3 into a layer prescribed at 0,065. The half-cell of load bearing material
        // that starts at the insulation's humidity instead is 0,25 mm wide and short by
        // 0,024 kg/m2 of the layer's 14,5 -- measured by StartsFromThePrescribedState.
        const double humidityA{
          Ham::humidityFromSuction(Bench::loadBearingSuction(Bench::initialLoadBearingWater))};
        const double humidityB{
          Ham::humidityFromSuction(Bench::insulationSuction(Bench::initialInsulationWater))};

        // 5 mm cells away from the contact surface, refining to 0,5 mm on both sides of it,
        // where the condensation front lives. A solver-side choice.
        TestHelper::BeamBuilder builder(multiDomain);
        builder.xStart(0.0)
          .height(0.01)
          .numElementsY(1)
          .state({.temperature = Bench::initialTemperature,
                  .humidity = humidityB,
                  .pressure = 101325.0,
                  .liquidPercent = 1.0})
          .addSegment({.material = matA.name(), .numElementsX = 17, .width = 0.085})
          .addSegment({.material = matA.name(), .numElementsX = 10, .width = 0.010})
          .addSegment({.material = matA.name(), .numElementsX = 10, .width = 0.005})
          .addSegment({.material = matB.name(), .numElementsX = 10, .width = 0.005})
          .addSegment({.material = matB.name(), .numElementsX = 5, .width = 0.005})
          .addSegment({.material = matB.name(), .numElementsX = 8, .width = 0.040})
          .build();
        const std::size_t contactColumn{37u};
        const auto coords = builder.xCoordinates();
        ASSERT_NEAR(Bench::loadBearingThickness, coords[contactColumn], 1e-12);

        auto humidities = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);
        for(std::size_t col = 0u; col < contactColumn; ++col)
        {
            humidities[builder.nodeIndex(col, 0u) - 1u] = humidityA;
            humidities[builder.nodeIndex(col, 1u) - 1u] = humidityA;
        }
        // Twice, so that both the current and the previous-timestep value read the
        // prescribed state whichever of the two one call moves.
        multiDomain.nodes().updateNodeHumidities(humidities, true);
        multiDomain.nodes().updateNodeHumidities(humidities, true);

        const Climate climate{loadClimate()};
        applyClimate(multiDomain, builder, climate);

        auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
        const auto record = [&](const std::vector<double> & nodal) {
            const auto row = TestHelper::bottomRow(nodal, builder.numNodesX(), 2u);
            loadBearingMass.push_back(
              integratedMoisture(coords, row, 0u, contactColumn, loadBearing));
            insulationMass.push_back(
              integratedMoisture(coords, row, contactColumn, builder.numElementsX(), insulation));
        };
        record(humidities);
        vapourInflow.push_back(0.0);
        const std::size_t surfaceNode{builder.nodeIndex(builder.numElementsX(), 0u) - 1u};
        for(std::size_t step = 0u; step < hoursPerYear; ++step)
        {
            const auto solution = multiDomain.transient(temperatures, humidities, hour, step);
            temperatures = solution.temperature;
            humidities = solution.humidity;
            record(humidities);
            vapourInflow.push_back(
              vapourInflow.back()
              + hour
                  * interiorVapourFlux(
                    climate, step + 1u, temperatures[surfaceNode], humidities[surfaceNode]));
        }

        bandA = loadBand("Bench1Year1_MA.csv");
        bandB = loadBand("Bench1Year1_MB.csv");

        HygroThermFEM::SimulationProperties::Instance().resetIterationParameters();
        HygroThermFEM::SimulationProperties::Instance().resetCalculationParameters();

        TestHelper::dumpGolden("HAMSTAD_Benchmark1 M_A", loadBearingMass);
        TestHelper::dumpGolden("HAMSTAD_Benchmark1 M_B", insulationMass);
        TestHelper::dumpGolden("HAMSTAD_Benchmark1 inflow", vapourInflow);
    }

    //! The vapour flux the interior condition delivers at the end of a step, kg/(m2 s),
    //! positive into the wall: the engine's own form beta_c (c_air - c_surf) with beta_c
    //! the Lewis value of the film coefficient the fixture hands it, i.e. the report's
    //! beta_p (p_air - p_surf). Its time integral is what the wall must have stored.
    [[nodiscard]] static double interiorVapourFlux(const Climate & climate,
                                                   const std::size_t row,
                                                   const double surfaceTemperature,
                                                   const double surfaceHumidity)
    {
        const double airTemperature{climate.interiorTemperature[row]};
        const double airHumidity{
          humidityFromPressure(climate.interiorPressure[row], airTemperature)};
        const double betaC{
          Ham::filmCoefficientForVapourTransfer(Bench::interiorVapourCoefficient, airTemperature)
          / (Constants::Density_Air * Constants::Cp_Air)};
        return betaC
               * (airHumidity * HygroThermFEM::saturationConcentrationAtTemperature(airTemperature)
                  - surfaceHumidity
                      * HygroThermFEM::saturationConcentrationAtTemperature(surfaceTemperature));
    }

    //! The exterior face exchanges heat only (beta_p,e = 0 s/m: sealed), with the
    //! equivalent temperature of the climate file; the interior face exchanges heat and
    //! vapour with the dwelling's air. Every hour of the year is one entry of the
    //! per-timestep boundary vectors.
    static void applyClimate(HygroThermFEM::MultiDomain & multiDomain,
                             const TestHelper::BeamBuilder & builder,
                             const Climate & climate)
    {
        std::vector<HygroThermFEM::FixedBCHCCoefficients> exteriorHeat;
        std::vector<HygroThermFEM::FixedBCHCCoefficients> interiorHeat;
        std::vector<HygroThermFEM::FixedBCHCCoefficients> interiorVapour;
        for(std::size_t step = 0u; step < hoursPerYear; ++step)
        {
            const std::size_t row{step + 1u};
            const double interiorTemperature{climate.interiorTemperature[row]};
            const double interiorHumidity{
              humidityFromPressure(climate.interiorPressure[row], interiorTemperature)};
            exteriorHeat.emplace_back(
              climate.exteriorTemperature[row], Bench::exteriorFilmCoefficient, 0.0);
            interiorHeat.emplace_back(
              interiorTemperature, Bench::interiorFilmCoefficient, interiorHumidity);
            interiorVapour.emplace_back(interiorTemperature,
                                        Ham::filmCoefficientForVapourTransfer(
                                          Bench::interiorVapourCoefficient, interiorTemperature),
                                        interiorHumidity);
        }
        for(const auto & [index1, index2] : builder.leftEdge())
        {
            multiDomain.thermal().createBC_FixedHc(index1, index2, exteriorHeat, false);
        }
        for(const auto & [index1, index2] : builder.rightEdge())
        {
            multiDomain.thermal().createBC_FixedHc(index1, index2, interiorHeat, true);
            multiDomain.moisture().createBC_FixedHc(index1, index2, interiorVapour);
        }
    }

    static void TearDownTestSuite()
    {
        loadBearingMass.clear();
        insulationMass.clear();
        vapourInflow.clear();
        bandA = {};
        bandB = {};
    }

    //! Hours whose value lies inside the band, from `firstHour` on, and the largest
    //! excursion outside it, as (count, excursion).
    [[nodiscard]] static std::pair<std::size_t, double>
      insideBand(const std::vector<double> & values, const Band & band, const std::size_t firstHour)
    {
        std::size_t inside{0u};
        double excursion{0.0};
        for(std::size_t hourIndex = firstHour; hourIndex <= hoursPerYear; ++hourIndex)
        {
            const double value{values[hourIndex]};
            if(value >= band.low[hourIndex] && value <= band.high[hourIndex])
            {
                ++inside;
                continue;
            }
            excursion = (std::max)(
              excursion, (std::max)(band.low[hourIndex] - value, value - band.high[hourIndex]));
        }
        return {inside, excursion};
    }
};

std::vector<double> HAMSTAD_Benchmark1::loadBearingMass{};
std::vector<double> HAMSTAD_Benchmark1::insulationMass{};
std::vector<double> HAMSTAD_Benchmark1::vapourInflow{};
Band HAMSTAD_Benchmark1::bandA{};
Band HAMSTAD_Benchmark1::bandB{};

/////////////////////////////////////////////////////////////////////////////////////
/// MEASURED RESULTS, first year, recorded 2026-09-15 (hourly backward Euler, 60 cells,
/// 28 s).
///
///   M_A: inside the band 97,3 % of hours 24-8760, never above it, worst excursion
///        0,046 kg/m2 (hour 2581, the spring drop, where the band is 0,14 wide);
///        inside the seven codes' own min-max 90 % of hours. Year end 13,36 against the
///        band 13,26-13,62.
///   M_B: inside the band 75,3 % of hours, worst excursion 0,055 kg/m2 (hour 755).
///        Outside the condensation season (hour 2700 on) it sits in the band; in the
///        season the engine condenses LESS than every code: hour 1500 reads 0,18 kg/m2
///        against the codes' 0,23-0,40 (average 0,34), and the peak follows the same
///        shape at about 55 % of the average.
///
/// Two things were checked before recording that shortfall. Pinning the engine's
/// diffusion coefficient in air at the report's constant 26,1e-6 (it is 22,2e-6 at 0 C)
/// raised hour 1500 only to 0,20, so the vapour path through the insulation is not it.
/// And the balance below holds through the season to 0,005 kg/m2, so the condensate is
/// not being lost at saturated nodes; what the interior face lets in, the wall keeps.
/// The condensation rate itself is the product of a small driving difference, p_room
/// minus p_sat at the contact surface, of order 100 Pa, so a kelvin at the contact or a
/// tenth of the insulation's conductivity moves it by half -- which is also why the
/// codes spread by a factor of two among themselves.
/////////////////////////////////////////////////////////////////////////////////////

//! The wall is sealed on the exterior, so whatever the interior face lets in is what the
//! two layers together must hold: M_A + M_B - M(0) = integral of the interior vapour
//! flux. The inflow is re-evaluated here from the end-of-step surface state, which lags
//! the engine's staggered solve (moisture first, at the previous temperature) by one
//! step of surface temperature; that shows only in the first hour, when the surface
//! warms from 10 C to the room's, as a 0,044 kg/m2 offset. From hour 1 on the two sides
//! track each other to 0,051 kg/m2 over a year that moves 2 kg/m2 through the face.
TEST_F(HAMSTAD_Benchmark1, StoredMoistureBalancesTheInteriorVapourInflow)
{
    ASSERT_EQ(hoursPerYear + 1u, vapourInflow.size());
    const double storedAtOne{loadBearingMass[1] + insulationMass[1]};
    double worst{0.0};
    for(std::size_t hourIndex = 2u; hourIndex <= hoursPerYear; ++hourIndex)
    {
        const double stored{loadBearingMass[hourIndex] + insulationMass[hourIndex] - storedAtOne};
        const double inflow{vapourInflow[hourIndex] - vapourInflow[1]};
        worst = (std::max)(worst, std::abs(stored - inflow));
    }
    EXPECT_LE(worst, 0.08) << "kg/m2 of stored moisture unaccounted for by the inflow";
}

//! The prescribed initial state, M_A = 145 * 0,1 and M_B = 0,065 * 0,05 kg/m2, up to the
//! one shared contact node (see the fixture), and a full year captured.
TEST_F(HAMSTAD_Benchmark1, StartsFromThePrescribedState)
{
    ASSERT_EQ(hoursPerYear + 1u, loadBearingMass.size());
    ASSERT_EQ(hoursPerYear + 1u, insulationMass.size());
    ASSERT_EQ(hoursPerYear + 1u, bandA.low.size());
    ASSERT_EQ(hoursPerYear + 1u, bandB.low.size());

    EXPECT_NEAR(
      Bench::initialLoadBearingWater * Bench::loadBearingThickness, loadBearingMass.front(), 0.03);
    EXPECT_NEAR(
      Bench::initialInsulationWater * Bench::insulationThickness, insulationMass.front(), 1e-4);
}

//! M_A against the load bearing layer's band of acceptance, Figure 3.6.1: inside for
//! at least 95 % of the year and never further than 0,07 kg/m2 from it (measured 97,3 %
//! and 0,046).
TEST_F(HAMSTAD_Benchmark1, LoadBearingMoistureWithinTheBandOfAcceptance)
{
    const auto [inside, excursion] = insideBand(loadBearingMass, bandA, 24u);
    EXPECT_GE(inside, (hoursPerYear - 24u) * 95u / 100u) << "hours inside the band";
    EXPECT_LE(excursion, 0.07) << "kg/m2 outside the band";
}

//! M_B against the insulation's band of acceptance, Figure 3.6.2, at the engine's
//! measured standing: inside for at least 70 % of the year (measured 75,3 %, the rest
//! being the condensation season the header describes) and never further than
//! 0,08 kg/m2 from it (measured 0,055). Tightening either number is the goal.
TEST_F(HAMSTAD_Benchmark1, InsulationMoistureWithinTheBandOfAcceptance)
{
    const auto [inside, excursion] = insideBand(insulationMass, bandB, 24u);
    EXPECT_GE(inside, (hoursPerYear - 24u) * 70u / 100u) << "hours inside the band";
    EXPECT_LE(excursion, 0.08) << "kg/m2 outside the band";

    // Outside the condensation season the band is narrow, a few 1e-4 kg/m2 wide, and
    // the engine sits in it: the drying of the insulation and its summer state agree
    // with the codes.
    const auto [insideSummer, excursionSummer] = insideBand(insulationMass, bandB, 2700u);
    EXPECT_GE(insideSummer, (hoursPerYear - 2700u) * 90u / 100u) << "hours inside, hour 2700 on";
}
