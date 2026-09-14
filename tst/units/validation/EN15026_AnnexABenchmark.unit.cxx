#include <algorithm>
#include <cstddef>
#include <map>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "EN15026Material.hxx"
#include "HygroThermFEM2D.hxx"
#include "SlabCreator.hxx"
#include "TestHelpers.hxx"

/////////////////////////////////////////////////////////////////////////////////////
/// EN 15026:2007 Annex A -- THE NORMATIVE BENCHMARK, against the standard's own
/// acceptance tables.
///
/// A.3 states the criterion: "The calculated results shall be within the limits given
/// in Tables A.1 and A.2." Those limits are transcribed verbatim below and are the
/// only expectations in this file. Nothing here is compared against a stored profile
/// of our own making.
///
/// THE RESULT: both tables satisfied in full, 19 of 19 water-content cells and 24 of
/// 24 temperature cells, on all three days.
///
/// Nearly every value lands within about 0,1 kg/m3 of its band CENTRE, against bands
/// 4,3 wide, so this is close agreement with the annex's analytic solution rather than
/// a scrape past its limits. The measured worst margin to any limit is 0,221 kg/m3 for
/// Table A.1 and 0,122 C for Table A.2.
///
/// WHAT MADE IT PASS, recorded because the failure is instructive. The annex's
/// diffusion resistance factor is moisture dependent, mu(w) = 200/f(w), running 212 at
/// the initial state to 866 at the wetted surface. The engine used to carry a single
/// scalar per material, so this benchmark ran at the annex's dry-range 200, resisted
/// vapour too little wherever the material was wet, and came out 5 of 19 cells too wet
/// by up to 2,0 kg/m3. That was a MODEL limit rather than a numerical one: doubling the
/// mesh and halving the time step moved the worst miss by 0,03 kg/m3. And no single
/// scalar could have fixed it -- at mu = 300 the profile was 16 of 19 and still
/// uniformly too wet, while at mu = 500 it fell outside in BOTH directions at once, too
/// wet at x = 0,02 and too dry at x = 0,06. The shape was wrong, not the level.
/// Supplying the annex's actual mu(w) curve is what closed it.
///
/// ONE ODDITY IN THE STANDARD, which this run appears to explain. Every printed band in
/// Table A.1 is 4,3 kg/m3 wide except day 365 at x = 0,04, printed 75,6 to 77,9, a width
/// of 2,3. Our value there is 75,821. Read with the printed minimum it sits 0,93 below
/// the band centre, far outside the roughly 0,1 spread of every other cell; read as if
/// the minimum were 73,6 -- restoring the 4,3 width and leaving the printed maximum
/// untouched -- it sits 0,07 from centre, exactly in line with the rest. A transposed
/// digit in the published table is the most economical explanation. The table is
/// transcribed verbatim below regardless, and the run passes either way.
/////////////////////////////////////////////////////////////////////////////////////

namespace
{
    namespace Annex = TestHelper::EN15026;

    //! One printed cell pair: the standard's validity limits at a day and a depth.
    struct Limit
    {
        int days;
        double position;   //!< [m]
        double minimum;
        double maximum;
    };

    //! Table A.1 (page 20) -- limits of validity for the humidity calculations,
    //! water content in kg/m3. Transcribed verbatim, blank cells omitted.
    //!
    //! Note day 365 at x = 0,04: the printed band is 75,6 to 77,9, a width of 2,3
    //! where every other cell in the table is 4,3 wide. That is what the standard
    //! prints and it is reproduced here unchanged.
    const std::vector<Limit> tableA1{
      {7, 0.01, 50.2, 54.5},     {7, 0.02, 41.3, 45.6},   {7, 0.03, 40.8, 45.1},
      {7, 0.04, 40.8, 45.1},     {30, 0.01, 81.0, 85.3},  {30, 0.02, 51.1, 55.3},
      {30, 0.03, 43.6, 47.9},    {30, 0.04, 41.5, 45.7},  {30, 0.05, 40.9, 45.2},
      {30, 0.06, 40.8, 45.1},    {30, 0.08, 40.8, 45.1},  {365, 0.01, 117.5, 121.8},
      {365, 0.02, 104.4, 108.7}, {365, 0.03, 88.7, 93.0}, {365, 0.04, 75.6, 77.9},
      {365, 0.05, 62.8, 67.1},   {365, 0.06, 55.7, 60.0}, {365, 0.08, 47.9, 52.2},
      {365, 0.10, 44.1, 48.4}};

    //! Table A.2 (page 21) -- limits of validity for the temperature calculations,
    //! degrees Celsius. Note the depths are METRES here while Table A.1's are
    //! centimetres: the thermal front runs far deeper than the moisture front, which
    //! is what makes the domain below tens of metres long.
    const std::vector<Limit> tableA2{
      {7, 0.5, 26.4, 26.9},   {7, 1.0, 23.6, 24.1},   {7, 1.5, 21.7, 22.2},
      {7, 2.0, 20.6, 21.1},   {7, 2.5, 20.0, 20.5},   {7, 3.0, 19.8, 20.4},
      {7, 4.0, 19.8, 20.3},   {7, 5.0, 19.8, 20.3},   {30, 0.5, 28.1, 28.6},
      {30, 1.0, 26.5, 27.0},  {30, 1.5, 25.0, 25.5},  {30, 2.0, 23.7, 24.3},
      {30, 2.5, 22.7, 23.2},  {30, 3.0, 21.8, 22.3},  {30, 4.0, 20.7, 21.2},
      {30, 5.0, 20.1, 20.6},  {365, 0.5, 29.2, 29.8}, {365, 1.0, 28.8, 29.3},
      {365, 1.5, 28.3, 28.8}, {365, 2.0, 27.8, 28.4}, {365, 2.5, 27.4, 27.9},
      {365, 3.0, 26.9, 27.4}, {365, 4.0, 26.0, 26.6}, {365, 5.0, 25.2, 25.7}};

    double profileAt(const std::vector<double> & coords,
                     const std::vector<double> & values,
                     const double position)
    {
        const auto upper = std::ranges::upper_bound(coords, position) - coords.begin();
        const auto right = (std::min)(static_cast<std::size_t>(upper), coords.size() - 1);
        const auto left = right - 1;
        const double fraction = (position - coords[left]) / (coords[right] - coords[left]);
        return values[left] + fraction * (values[right] - values[left]);
    }
}   // namespace

//! Runs the benchmark ONCE for the whole suite: it is the most expensive test in the
//! repository and both enabled tests read the same three profiles from it.
class EN15026_AnnexABenchmark : public ::testing::Test
{
protected:
    static std::vector<double> coordinates;
    static std::map<int, std::vector<double>> waterByDay;
    static std::map<int, std::vector<double>> temperatureByDay;

    static void SetUpTestSuite()
    {
        // Every physical term in, plus the moisture-dependent conductivity branch:
        // the annex specifies lambda(w), liquid transport and latent heat together.
        HygroThermFEM::SimulationProperties::Instance().setCalculationParameters(
          false, false, false, false, true);
        HygroThermFEM::SimulationProperties::Instance().setIterationParameters(1.0, 1e-9, 500);

        HygroThermFEM::MultiDomain multiDomain({.performThermal = true, .performMoisture = true});

        // The material the annex actually specifies, resistance factor curve included.
        const auto & material = multiDomain.materials().createSolidMaterial(
          Annex::material(Annex::dryResistanceFactor, true));

        // "Semi-infinite" in practice: Table A.2 reads temperature out to 5 m, and over
        // 365 days the thermal front reaches roughly 12 m, so the far end must sit well
        // beyond that to stay undisturbed. Cells grow geometrically from half a
        // millimetre at the wetted surface, where the moisture gradient is steep.
        coordinates = Annex::gradedCoordinates(5.0e-4, 1.06, 30.0);
        const auto nColumns = coordinates.size();

        TestHelper::SlabBuilder(multiDomain)
          .gridXCoordinates(coordinates)
          .height(0.05)
          .material(material.name())
          .state({.temperature = 20.0,
                  .humidity = 0.5,
                  .pressure = 101325.0,
                  .liquidPercent = 1.0})
          .build();

        // A.2: uniform 20 C and phi = 0,5, then the surface steps to 30 C and phi = 0,95
        // with no boundary resistance, so both fields are imposed directly on the edge.
        multiDomain.thermal().createBC_FixedTemperature(1, 2, 30.0);
        multiDomain.moisture().createBC_FixedHumidity(
          1, 2, HygroThermFEM::TemperatureAndHumidity{30.0, 0.95});

        // Backward Euler, refined through the first hour where the step change is
        // sharpest, then six-hourly to a year. Block boundaries land exactly on the
        // three days the annex asks for. Six-hourly is a cost choice, not an accuracy
        // one: two-hourly steps move the worst Table A.1 miss by 0,007 kg/m3 and
        // daily steps by 0,03, against bands 4,3 wide.
        const std::vector<std::pair<double, int>> timeBlocks{
          {30.0, 120}, {600.0, 138}, {21600.0, 24}, {21600.0, 92}, {21600.0, 1340}};
        const std::map<int, int> captureAt{{282, 7}, {374, 30}, {1714, 365}};

        auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
        auto humidities = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);

        int step{0};
        for(const auto & [dTime, nSteps] : timeBlocks)
        {
            for(int blockStep = 0; blockStep < nSteps; ++blockStep)
            {
                const auto solution = multiDomain.transient(temperatures, humidities, dTime, step);
                temperatures = solution.temperature;
                humidities = solution.humidity;
                ++step;
                const auto capture = captureAt.find(step);
                if(capture != captureAt.end())
                {
                    waterByDay[capture->second] =
                      TestHelper::bottomRow(solution.waterContent, nColumns, 2);
                    temperatureByDay[capture->second] =
                      TestHelper::bottomRow(temperatures, nColumns, 2);
                }
            }
        }

        HygroThermFEM::SimulationProperties::Instance().resetIterationParameters();
        HygroThermFEM::SimulationProperties::Instance().resetCalculationParameters();
    }

    static void TearDownTestSuite()
    {
        coordinates.clear();
        waterByDay.clear();
        temperatureByDay.clear();
    }

    //! The computed value at one printed cell.
    static double valueAt(const std::map<int, std::vector<double>> & field, const Limit & limit)
    {
        return profileAt(coordinates, field.at(limit.days), limit.position);
    }
};

std::vector<double> EN15026_AnnexABenchmark::coordinates{};
std::map<int, std::vector<double>> EN15026_AnnexABenchmark::waterByDay{};
std::map<int, std::vector<double>> EN15026_AnnexABenchmark::temperatureByDay{};

//! All three days were captured, the surface pins held, and the far end stayed
//! undisturbed -- so the semi-infinite region the benchmark assumes is still a fair
//! description at a year. Checked first, because every limit below is meaningless if
//! the domain ran out.
TEST_F(EN15026_AnnexABenchmark, BenchmarkRanAsTheAnnexSpecifies)
{
    ASSERT_EQ(3u, waterByDay.size());
    ASSERT_EQ(3u, temperatureByDay.size());

    for(const int day : {7, 30, 365})
    {
        EXPECT_NEAR(30.0, temperatureByDay.at(day).front(), 1e-9) << "surface pin, day " << day;

        // A hundredth of a degree at 30 m by day 365 (measured 0,003). Table A.2 reads
        // no deeper than 5 m, where the profile is still 25 C, so a disturbance this
        // small at the outer boundary cannot reach the cells being judged.
        EXPECT_NEAR(20.0, temperatureByDay.at(day).back(), 1e-2)
          << "the domain is too short: the thermal front reached the far end by day " << day;
    }

    // The moisture front is orders of magnitude slower than the thermal one, so the far
    // end must not move AT ALL. Compared day against day rather than against the
    // isotherm, because the engine reads water content from a table: the far end sits
    // at 42,981 where the closed form gives 42,922, and that 0,06 kg/m3 is the
    // tabulation offset measured in EN15026_AnnexAProperties, present from step one.
    EXPECT_NEAR(waterByDay.at(7).back(), waterByDay.at(365).back(), 1e-6)
      << "the moisture front reached the far end between day 7 and day 365";
    EXPECT_NEAR(Annex::waterContent(0.5), waterByDay.at(365).back(), 0.1)
      << "the undisturbed far end is not at the benchmark's initial state";
}

//! EN 15026:2007 A.3 against Table A.2. This one PASSES outright: all 24 printed
//! cells, on all three days, at depths from 0,5 m to 5 m.
TEST_F(EN15026_AnnexABenchmark, TemperatureProfilesSatisfyTableA2)
{
    for(const auto & limit : tableA2)
    {
        const double value{valueAt(temperatureByDay, limit)};
        EXPECT_GE(value, limit.minimum)
          << "day " << limit.days << ", x = " << limit.position << " m: " << value
          << " C is below the standard's minimum " << limit.minimum;
        EXPECT_LE(value, limit.maximum)
          << "day " << limit.days << ", x = " << limit.position << " m: " << value
          << " C is above the standard's maximum " << limit.maximum;
    }
}

/////////////////////////////////////////////////////////////////////////////////////
/// MEASURED RESULTS, cell by cell, recorded 2026-09-14.
///
/// "Centre" is the midpoint of the printed band, which is the annex's analytic solution
/// -- the limits are that solution plus and minus 2,5 % of the profile's range. So the
/// last column is the distance from the standard's own answer, and it is the number
/// worth watching: the pass/fail margin only says we are inside a fairly wide band.
///
///   Table A.1 -- water content [kg/m3]
///    day     x       ours       band        centre    diff    result
///      7  0,01     52,611   50,2 - 54,5     52,35    +0,26    PASS
///      7  0,02     43,544   41,3 - 45,6     43,45    +0,09    PASS
///      7  0,03     43,001   40,8 - 45,1     42,95    +0,05    PASS
///      7  0,04     42,983   40,8 - 45,1     42,95    +0,03    PASS
///     30  0,01     83,282   81,0 - 85,3     83,15    +0,13    PASS
///     30  0,02     53,219   51,1 - 55,3     53,20    +0,02    PASS
///     30  0,03     45,777   43,6 - 47,9     45,75    +0,03    PASS
///     30  0,04     43,622   41,5 - 45,7     43,60    +0,02    PASS
///     30  0,05     43,111   40,9 - 45,2     43,05    +0,06    PASS
///     30  0,06     43,002   40,8 - 45,1     42,95    +0,05    PASS
///     30  0,08     42,983   40,8 - 45,1     42,95    +0,03    PASS
///    365  0,01    119,749  117,5 - 121,8   119,65    +0,10    PASS
///    365  0,02    106,660  104,4 - 108,7   106,55    +0,11    PASS
///    365  0,03     90,942   88,7 - 93,0     90,85    +0,09    PASS
///    365  0,04     75,821   75,6 - 77,9     76,75    -0,93    PASS  (narrow band)
///    365  0,05     64,891   62,8 - 67,1     64,95    -0,06    PASS
///    365  0,06     57,742   55,7 - 60,0     57,85    -0,11    PASS
///    365  0,08     49,998   47,9 - 52,2     50,05    -0,05    PASS
///    365  0,10     46,273   44,1 - 48,4     46,25    +0,02    PASS
///                                                  19 of 19 PASS
///
///   Table A.2 -- temperature [C]
///    day     x       ours       band        centre    diff    result
///      7  0,5      26,555   26,4 - 26,9     26,65    -0,10    PASS
///      7  1,0      23,723   23,6 - 24,1     23,85    -0,13    PASS
///      7  1,5      21,822   21,7 - 22,2     21,95    -0,13    PASS
///      7  2,0      20,767   20,6 - 21,1     20,85    -0,08    PASS
///      7  2,5      20,279   20,0 - 20,5     20,25    +0,03    PASS
///      7  3,0      20,089   19,8 - 20,4     20,10    -0,01    PASS
///      7  4,0      20,006   19,8 - 20,3     20,05    -0,04    PASS
///      7  5,0      20,000   19,8 - 20,3     20,05    -0,05    PASS
///     30  0,5      28,316   28,1 - 28,6     28,35    -0,03    PASS
///     30  1,0      26,690   26,5 - 27,0     26,75    -0,06    PASS
///     30  1,5      25,207   25,0 - 25,5     25,25    -0,04    PASS
///     30  2,0      23,916   23,7 - 24,3     24,00    -0,08    PASS
///     30  2,5      22,841   22,7 - 23,2     22,95    -0,11    PASS
///     30  3,0      21,988   21,8 - 22,3     22,05    -0,06    PASS
///     30  4,0      20,869   20,7 - 21,2     20,95    -0,08    PASS
///     30  5,0      20,328   20,1 - 20,6     20,35    -0,02    PASS
///    365  0,5      29,528   29,2 - 29,8     29,50    +0,03    PASS
///    365  1,0      29,040   28,8 - 29,3     29,05    -0,01    PASS
///    365  1,5      28,555   28,3 - 28,8     28,55    +0,01    PASS
///    365  2,0      28,076   27,8 - 28,4     28,10    -0,02    PASS
///    365  2,5      27,604   27,4 - 27,9     27,65    -0,05    PASS
///    365  3,0      27,140   26,9 - 27,4     27,15    -0,01    PASS
///    365  4,0      26,246   26,0 - 26,6     26,30    -0,05    PASS
///    365  5,0      25,404   25,2 - 25,7     25,45    -0,05    PASS
///                                                  24 of 24 PASS
///
///   Worst |diff|: 0,26 kg/m3 on Table A.1 against bands 4,3 wide, and 0,13 C on
///   Table A.2 against bands 0,5 wide. Every Table A.1 entry other than the narrow-band
///   cell sits within 0,26 of the annex's own answer.
///
/// The -0,93 at day 365, x = 0,04 is the cell whose printed band is half the usual
/// width; measured against a band restored to 4,3 it would read -0,07, in line with
/// every other entry. See the note at tableA1 below.
/////////////////////////////////////////////////////////////////////////////////////

//! EN 15026:2007 A.3 against Table A.1. Every printed cell, on all three days.
//!
//! Worst measured margin to a limit is 0,221 kg/m3, at the one cell whose printed band
//! is half the width of every other; the next tightest is 1,889. See this file's header.
TEST_F(EN15026_AnnexABenchmark, WaterContentProfilesSatisfyTableA1)
{
    for(const auto & limit : tableA1)
    {
        const double value{valueAt(waterByDay, limit)};
        EXPECT_GE(value, limit.minimum)
          << "day " << limit.days << ", x = " << limit.position << " m: " << value
          << " kg/m3 is below the standard's minimum " << limit.minimum;
        EXPECT_LE(value, limit.maximum)
          << "day " << limit.days << ", x = " << limit.position << " m: " << value
          << " kg/m3 is above the standard's maximum " << limit.maximum;
    }
}
