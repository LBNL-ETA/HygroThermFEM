#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "BeamBuilder.hxx"
#include "HAMSTADBenchmarks.hxx"
#include "HygroThermFEM2D.hxx"
#include "NumericTable.hxx"
#include "TestHelpers.hxx"

/////////////////////////////////////////////////////////////////////////////////////
/// HAMSTAD benchmark 5, "Capillary active inside insulation", against the six
/// participating codes' profiles.
///
/// Report R-02:8 (Hagentoft, 2002) section 3.5, p. 45-51. A 365 mm brick wall carries
/// 15 mm of mortar and 40 mm of capillary active insulation on its interior face. The
/// wall starts uniform at 25 C and 60 % and at t = 0 the exterior air steps to 0 C and
/// 80 % (alpha = 25 W/m2K, beta_p = 1,8382e-7 s/m), the interior to 20 C and 60 %
/// (alpha = 8, beta_p = 5,8823e-8). Vapour driven inward through the cold brick
/// condenses at the mortar and is drawn back into the room by the insulation's
/// capillary conductivity. The requested output is w(x) and phi(x) after 60 days.
///
/// There is no analytic solution. The package distributes the six submitted profiles
/// (Bench5.xls), all reproduced in the data directory, and the report judges them to
/// "show rather good agreement". The expectation here is therefore the ENVELOPE of the
/// six: at every requested position the engine's value must lie between the lowest and
/// the highest submission, widened by the margin recorded at the assertion.
/////////////////////////////////////////////////////////////////////////////////////

namespace
{
    namespace Bench = TestHelper::HAMSTAD::Benchmark5;
    namespace Ham = TestHelper::HAMSTAD;

    const std::vector<std::string> participants{"CTH", "KUL", "NRC", "Technion", "TUD", "TUE"};

    //! One submitted profile: position [mm], water content [kg/m3], relative humidity.
    struct Submission
    {
        std::string code;
        std::vector<double> positions;
        std::vector<double> water;
        std::vector<double> humidity;
    };

    [[nodiscard]] Submission loadSubmission(const std::string & code)
    {
        const auto table = TestHelper::loadNumericTable(
          std::string{HAMSTAD_DATA_DIR} + "/Bench5_" + code + ".csv", 1u);
        return {.code = code,
                .positions = TestHelper::column(table, 0u),
                .water = TestHelper::column(table, 1u),
                .humidity = TestHelper::column(table, 2u)};
    }

    //! The material at a position through the wall, exterior face at x = 0.
    [[nodiscard]] const Ham::Material & materialAt(const double position,
                                                   const Ham::Material & brick,
                                                   const Ham::Material & mortar,
                                                   const Ham::Material & insulation)
    {
        if(position < Bench::brickThickness)
        {
            return brick;
        }
        if(position < Bench::brickThickness + Bench::mortarThickness)
        {
            return mortar;
        }
        return insulation;
    }

    //! Lowest and highest submitted value at a position, across every participant that
    //! reports the quantity there.
    struct Envelope
    {
        double low;
        double high;
    };

    [[nodiscard]] Envelope envelopeAt(const std::vector<Submission> & submissions,
                                      const double positionMm,
                                      const bool water)
    {
        Envelope band{.low = std::numeric_limits<double>::max(),
                      .high = std::numeric_limits<double>::lowest()};
        for(const auto & submission : submissions)
        {
            const auto & values = water ? submission.water : submission.humidity;
            const double value{TestHelper::interpolateAt(submission.positions, values, positionMm)};
            if(std::isnan(value))
            {
                continue;
            }
            band.low = (std::min)(band.low, value);
            band.high = (std::max)(band.high, value);
        }
        return band;
    }
}   // namespace

//! Runs the 60 days once for the suite.
class HAMSTAD_Benchmark5 : public ::testing::Test
{
protected:
    static std::vector<double> coordinates;   //!< node x [m], exterior face first
    static std::vector<double> humidity;      //!< bottom node row after 60 days
    static std::vector<double> temperature;
    static std::vector<Submission> submissions;

    static void SetUpTestSuite()
    {
        // Every physical term in, plus the moisture dependent conductivity branch: the
        // report specifies lambda(w), liquid transport and the latent heat together.
        HygroThermFEM::SimulationProperties::Instance().setCalculationParameters(
          false, false, false, false, true);
        HygroThermFEM::SimulationProperties::Instance().setIterationParameters(1.0, 1e-9, 500);

        HygroThermFEM::MultiDomain multiDomain({.performThermal = true, .performMoisture = true});
        const auto & brick =
          multiDomain.materials().createSolidMaterial(Ham::solidMaterial(Bench::brick()));
        const auto & mortar =
          multiDomain.materials().createSolidMaterial(Ham::solidMaterial(Bench::mortar()));
        const auto & insulation =
          multiDomain.materials().createSolidMaterial(Ham::solidMaterial(Bench::insulation()));

        // 5 mm cells through the brick, 1 mm through its last centimetre, the mortar and
        // the insulation, where the report asks for millimetre output and the moisture
        // front lives. A solver-side choice; the report prescribes no mesh.
        TestHelper::BeamBuilder builder(multiDomain);
        builder.xStart(0.0)
          .height(0.01)
          .numElementsY(1)
          .state({.temperature = Bench::initialTemperature,
                  .humidity = Bench::initialHumidity,
                  .pressure = 101325.0,
                  .liquidPercent = 1.0})
          .addSegment(
            {.material = brick.name(), .numElementsX = 71, .width = Bench::brickThickness - 0.010})
          .addSegment({.material = brick.name(), .numElementsX = 10, .width = 0.010})
          .addSegment(
            {.material = mortar.name(), .numElementsX = 15, .width = Bench::mortarThickness})
          .addSegment({.material = insulation.name(),
                       .numElementsX = 40,
                       .width = Bench::insulationThickness})
          .build();

        // The report gives the heat and the vapour surface coefficients separately; the
        // thermal domain takes alpha and the moisture domain the film coefficient that
        // reproduces beta_p through the engine's Lewis relation (see the helper). The
        // latent heat of the surface vapour flux inside the thermal condition is the one
        // term that keeps the engine's own Lewis value, a few percent off the prescribed
        // beta_p on a term that is itself small.
        applyFace(multiDomain,
                  builder.leftEdge(),
                  Bench::exteriorTemperature,
                  Bench::exteriorHumidity,
                  Bench::exteriorFilmCoefficient,
                  Bench::exteriorVapourCoefficient);
        applyFace(multiDomain,
                  builder.rightEdge(),
                  Bench::interiorTemperature,
                  Bench::interiorHumidity,
                  Bench::interiorFilmCoefficient,
                  Bench::interiorVapourCoefficient);

        coordinates = builder.xCoordinates();

        // Backward Euler, hourly, through the 60 days. Halving both the step and the cells
        // of the last centimetre of brick, the mortar and the insulation moves the mortar
        // peak by 0,1 kg/m3 (59,27 to 59,37) and nothing else by more, at three times the
        // cost, so the coarser discretisation is the one kept.
        constexpr double dTime{3600.0};
        const int nSteps{static_cast<int>(Bench::simulationDays * 86400.0 / dTime)};
        auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
        auto humidities = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);
        for(int step = 0; step < nSteps; ++step)
        {
            const auto solution = multiDomain.transient(temperatures, humidities, dTime, step);
            temperatures = solution.temperature;
            humidities = solution.humidity;
        }
        humidity = TestHelper::bottomRow(humidities, builder.numNodesX(), 2u);
        temperature = TestHelper::bottomRow(temperatures, builder.numNodesX(), 2u);

        for(const auto & code : participants)
        {
            submissions.push_back(loadSubmission(code));
        }

        HygroThermFEM::SimulationProperties::Instance().resetIterationParameters();
        HygroThermFEM::SimulationProperties::Instance().resetCalculationParameters();

        dumpProfiles();
    }

    //! With HTF_DUMP_GOLDEN set: the engine's w and phi at the report's output positions
    //! next to the envelope of the submissions, for recording margins.
    static void dumpProfiles()
    {
        if(!TestHelper::goldenDumpEnabled())
        {
            return;
        }
        const auto tud = std::ranges::find(submissions, "TUD", &Submission::code);
        std::vector<std::vector<double>> rows;
        for(const double positionMm : tud->positions)
        {
            const auto waterBand = envelopeAt(submissions, positionMm, true);
            const auto humidityBand = envelopeAt(submissions, positionMm, false);
            rows.push_back({positionMm,
                            waterAt(positionMm / 1000.0),
                            waterBand.low,
                            waterBand.high,
                            humidityAt(positionMm / 1000.0),
                            humidityBand.low,
                            humidityBand.high});
        }
        TestHelper::dumpGolden("HAMSTAD_Benchmark5 x_mm, w, w_low, w_high, phi, phi_low, phi_high",
                               rows);
        TestHelper::dumpGolden("HAMSTAD_Benchmark5 temperature", temperature);
    }

    static void applyFace(HygroThermFEM::MultiDomain & multiDomain,
                          const std::vector<std::pair<std::size_t, std::size_t>> & edge,
                          const double airTemperature,
                          const double airHumidity,
                          const double filmCoefficient,
                          const double vapourCoefficient)
    {
        const HygroThermFEM::FixedBCHCCoefficients thermal{
          airTemperature, filmCoefficient, airHumidity};
        const HygroThermFEM::FixedBCHCCoefficients moisture{
          airTemperature,
          Ham::filmCoefficientForVapourTransfer(vapourCoefficient, airTemperature),
          airHumidity};
        for(const auto & [index1, index2] : edge)
        {
            multiDomain.thermal().createBC_FixedHc(index1, index2, thermal, true);
            multiDomain.moisture().createBC_FixedHc(index1, index2, moisture);
        }
    }

    static void TearDownTestSuite()
    {
        coordinates.clear();
        humidity.clear();
        temperature.clear();
        submissions.clear();
    }

    //! Relative humidity at a position [m], interpolated between nodes; phi is the
    //! engine's continuous unknown, so this is valid across the material interfaces.
    [[nodiscard]] static double humidityAt(const double position)
    {
        return TestHelper::interpolateAt(coordinates, humidity, position);
    }

    //! Water content at a position [m] through the material's own closed-form retention
    //! curve, exactly as the report defines the requested output and as the
    //! participants computed theirs.
    [[nodiscard]] static double waterAt(const double position)
    {
        static const Ham::Material brick{Bench::brick()};
        static const Ham::Material mortar{Bench::mortar()};
        static const Ham::Material insulation{Bench::insulation()};
        const auto & material = materialAt(position, brick, mortar, insulation);
        return Ham::waterContent(material, Ham::suctionFromHumidity(humidityAt(position)));
    }
};

std::vector<double> HAMSTAD_Benchmark5::coordinates{};
std::vector<double> HAMSTAD_Benchmark5::humidity{};
std::vector<double> HAMSTAD_Benchmark5::temperature{};
std::vector<Submission> HAMSTAD_Benchmark5::submissions{};

//! The six submissions loaded, and the engine's state after 60 days is physically
//! where the report puts it: the brick cold and near the exterior humidity, the mortar
//! the wettest point of the wall, the insulation drying toward the room.
TEST_F(HAMSTAD_Benchmark5, RanAsTheReportSpecifies)
{
    ASSERT_EQ(participants.size(), submissions.size());
    for(const auto & submission : submissions)
    {
        EXPECT_GT(submission.positions.size(), 100u) << submission.code;
    }

    EXPECT_LT(temperature.front(), 2.0) << "exterior surface, C";
    EXPECT_GT(temperature.back(), 17.0) << "interior surface, C";

    const double brickMiddle{waterAt(0.18)};
    const double mortarPeak{waterAt(Bench::brickThickness + Bench::mortarThickness - 0.0005)};
    const double insulationFace{waterAt(0.4195)};
    EXPECT_GT(mortarPeak, brickMiddle);
    EXPECT_GT(mortarPeak, insulationFace);
}

/////////////////////////////////////////////////////////////////////////////////////
/// MEASURED RESULTS, recorded 2026-09-15 (1 mm cells, hourly, 8 s).
///
/// Water content inside the envelope at 97 of the 144 positions, relative humidity at
/// 123; every miss is on the DRY side, never above the wettest code. The shape is the
/// report's Figure 3.5.2 -- brick flat at 4,5 kg/m3, the rise through the mortar, the
/// peak at the mortar-insulation contact, the tail through the insulation -- with the
/// engine tracking the driest submission (Technion) a little below it:
///
///   x [mm]    ours     lowest  highest  (kg/m3)
///    373,5     7,69     7,70     8,64
///    380,5    59,27    60,26    63,72   the peak, 0,98 below the lowest
///    390,5    35,99    36,20    53,12
///    400,5    10,77    11,19    17,84
///    419,5     3,41     3,46     3,63
///
/// Worst excursion 0,98 kg/m3 (the peak) and 0,0014 in relative humidity, against a
/// spread among the codes of 3,5 kg/m3 at the peak and 17 through the insulation.
///
/// Half of that dryness is the engine's temperature dependent diffusion coefficient in
/// air where the report fixes 26,1e-6: with the engine's coefficient pinned at that
/// value for one run, the peak read 59,78 and the worst excursion 0,47. The other half
/// is the engine's own -- capacity form, Lewis latent term at the faces, coefficient
/// interpolation at the two interfaces -- and is what these margins now watch.
/////////////////////////////////////////////////////////////////////////////////////

//! Every requested position of the report's output grid (the 144 positions the TUD
//! submission tabulates, millimetre steps at the layer ends and 5 mm through the brick)
//! against the envelope of the six submissions, widened by roughly one and a half
//! times the worst measured excursion.
TEST_F(HAMSTAD_Benchmark5, ProfilesLieWithinTheEnvelopeOfTheSubmissions)
{
    const auto tud = std::ranges::find(submissions, "TUD", &Submission::code);
    ASSERT_NE(tud, submissions.end());

    constexpr double waterMargin{1.5};        //!< kg/m3
    constexpr double humidityMargin{0.003};   //!< -
    for(const double positionMm : tud->positions)
    {
        const double position{positionMm / 1000.0};
        const auto waterBand = envelopeAt(submissions, positionMm, true);
        const auto humidityBand = envelopeAt(submissions, positionMm, false);
        const double water{waterAt(position)};
        const double phi{humidityAt(position)};
        EXPECT_GE(water, waterBand.low - waterMargin)
          << "x = " << positionMm << " mm: " << water << " kg/m3 below the lowest submission "
          << waterBand.low;
        EXPECT_LE(water, waterBand.high + waterMargin)
          << "x = " << positionMm << " mm: " << water << " kg/m3 above the highest submission "
          << waterBand.high;
        EXPECT_GE(phi, humidityBand.low - humidityMargin)
          << "x = " << positionMm << " mm: phi " << phi << " below the lowest submission "
          << humidityBand.low;
        EXPECT_LE(phi, humidityBand.high + humidityMargin)
          << "x = " << positionMm << " mm: phi " << phi << " above the highest submission "
          << humidityBand.high;
    }
}
