#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "HAMSTADBenchmarks.hxx"
#include "HygroThermFEM2D.hxx"
#include "SlabCreator.hxx"
#include "TestHelpers.hxx"

/////////////////////////////////////////////////////////////////////////////////////
/// HAMSTAD benchmark 2, "Analytical benchmark" -- the one case of the package with a
/// closed-form solution, against the package's own tabulation of it.
///
/// Report R-02:8 (Hagentoft, 2002) section 3.2, p. 19-23. A 0,2 m homogeneous layer,
/// isothermal at 20 C, starts in equilibrium with air at 95 % relative humidity
/// (w = 84,7687 kg/m3). At t = 0 the air drops to 45 % on the exterior face (x = 0) and
/// 65 % on the interior face (x = 0,2 m) behind a surface coefficient of 1e-3 s/m,
/// which is no resistance at all: the printed solution pins both faces at their
/// equilibrium water contents from the first instant. The material has a constant
/// moisture diffusivity of 6e-10 m2/s, so w(x, t) obeys the linear diffusion equation
/// and the requested output is w(x) at 100, 300 and 1000 hours, at 0,01 m intervals.
///
/// The package distributes the analytical values in Bench2.xls next to the seven
/// participating codes' results, and that column is the only expectation here.
/////////////////////////////////////////////////////////////////////////////////////

namespace
{
    namespace Bench = TestHelper::HAMSTAD::Benchmark2;

    //! Bench2.xls, sheet "Page1", column "Analytical": water content [kg/m3] at
    //! x = 0 to 0,20 m in 0,01 m steps after 100, 300 and 1000 hours, four decimals,
    //! transcribed verbatim.
    struct AnalyticalRow
    {
        double position;   //!< [m]
        double hour100;
        double hour300;
        double hour1000;
    };

    const std::vector<AnalyticalRow> analytical{
      {0.00, 19.5356, 19.5356, 19.5356}, {0.01, 43.6440, 33.8097, 27.2914},
      {0.02, 62.8553, 47.0304, 34.8419}, {0.03, 75.0545, 58.3716, 41.9918},
      {0.04, 81.2271, 67.3823, 48.5646}, {0.05, 83.7155, 74.0123, 54.4092},
      {0.06, 84.5148, 78.5282, 59.4036}, {0.07, 84.7193, 81.3703, 63.4561},
      {0.08, 84.7610, 83.0086, 66.5036}, {0.09, 84.7677, 83.8367, 68.5088},
      {0.10, 84.7685, 84.1147, 69.4563}, {0.11, 84.7679, 83.9483, 69.3502},
      {0.12, 84.7623, 83.2873, 68.2114}, {0.13, 84.7276, 81.9357, 66.0777},
      {0.14, 84.5575, 79.5757, 63.0043}, {0.15, 83.8926, 75.8206, 59.0659},
      {0.16, 81.8227, 70.3060, 54.3585}, {0.17, 76.6882, 62.8107, 49.0001},
      {0.18, 66.5405, 53.3768, 43.1305}, {0.19, 50.5599, 42.3795, 36.9083},
      {0.20, 30.5059, 30.5059, 30.5059}};

    //! Uniform 1 mm cells across the 0,2 m layer. A solver-side choice; the report
    //! prescribes no mesh.
    [[nodiscard]] std::vector<double> uniformCoordinates(const std::size_t cells,
                                                         const double length)
    {
        std::vector<double> coords;
        coords.reserve(cells + 1u);
        for(std::size_t index = 0u; index <= cells; ++index)
        {
            coords.push_back(length * static_cast<double>(index) / static_cast<double>(cells));
        }
        return coords;
    }

    //! The value at a printed position, read off the bottom node row. Every printed
    //! position lands on a node, so this is a lookup rather than an interpolation.
    [[nodiscard]] double valueAt(const std::vector<double> & coords,
                                 const std::vector<double> & values,
                                 const double position)
    {
        const auto nearest = std::ranges::min_element(
          coords, {}, [position](const double coord) { return std::abs(coord - position); });
        return values[static_cast<std::size_t>(nearest - coords.begin())];
    }
}   // namespace

//! Runs the benchmark once for the suite and keeps the three requested profiles.
class HAMSTAD_Benchmark2 : public ::testing::Test
{
protected:
    static std::vector<double> coordinates;
    static std::map<int, std::vector<double>> waterByHour;
    static std::vector<double> initialWater;

    static void SetUpTestSuite()
    {
        HygroThermFEM::SimulationProperties::Instance().setCalculationParameters(
          false, false, false, false, false);
        HygroThermFEM::SimulationProperties::Instance().setIterationParameters(1.0, 1e-9, 500);

        // Isothermal by definition ("decoupling of thermal and moisture process"), so only
        // the moisture domain is solved and every node stays at 20 C.
        HygroThermFEM::MultiDomain multiDomain({.performThermal = false, .performMoisture = true});
        const auto & material = multiDomain.materials().createSolidMaterial(Bench::material());

        coordinates = uniformCoordinates(200u, Bench::thickness);
        const auto nColumns = coordinates.size();
        TestHelper::SlabBuilder(multiDomain)
          .gridXCoordinates(coordinates)
          .height(0.01)
          .material(material.name())
          .state({.temperature = Bench::temperature,
                  .humidity = Bench::initialHumidity,
                  .pressure = 101325.0,
                  .liquidPercent = 1.0})
          .build();

        // The prescribed beta_p = 1e-3 s/m on both faces, handed to the engine as the film
        // coefficient that reproduces it through its Lewis relation (see the helper).
        const double film{TestHelper::HAMSTAD::filmCoefficientForVapourTransfer(
          Bench::vapourCoefficient, Bench::temperature)};
        const std::size_t lastColumn{2u * nColumns};
        multiDomain.moisture().createBC_FixedHc(
          1u,
          2u,
          HygroThermFEM::FixedBCHCCoefficients{Bench::temperature, film, Bench::exteriorHumidity});
        multiDomain.moisture().createBC_FixedHc(
          lastColumn - 1u,
          lastColumn,
          HygroThermFEM::FixedBCHCCoefficients{Bench::temperature, film, Bench::interiorHumidity});

        // Backward Euler: ten-minute steps through the first ten hours, where the faces
        // have just dropped, then hourly to 1000 h. Block boundaries land on the three
        // requested hours.
        const std::vector<std::pair<double, int>> timeBlocks{{600.0, 60}, {3600.0, 990}};
        const std::map<int, int> captureAt{{150, 100}, {350, 300}, {1050, 1000}};

        auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
        auto humidities = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);
        initialWater = TestHelper::bottomRow(
          multiDomain.nodes().properties(HygroThermFEM::Variable::water), nColumns, 2u);

        int step{0};
        for(const auto & [dTime, nSteps] : timeBlocks)
        {
            for(int blockStep = 0; blockStep < nSteps; ++blockStep)
            {
                const auto solution = multiDomain.transient(temperatures, humidities, dTime, step);
                temperatures = solution.temperature;
                humidities = solution.humidity;
                ++step;
                if(const auto capture = captureAt.find(step); capture != captureAt.end())
                {
                    waterByHour[capture->second] =
                      TestHelper::bottomRow(solution.waterContent, nColumns, 2u);
                }
            }
        }

        HygroThermFEM::SimulationProperties::Instance().resetIterationParameters();
        HygroThermFEM::SimulationProperties::Instance().resetCalculationParameters();

        for(const auto & [hour, water] : waterByHour)
        {
            std::vector<double> atPrintedPositions;
            for(const auto & row : analytical)
            {
                atPrintedPositions.push_back(valueAt(coordinates, water, row.position));
            }
            TestHelper::dumpGolden("HAMSTAD_Benchmark2 hour " + std::to_string(hour),
                                   atPrintedPositions);
        }
    }

    static void TearDownTestSuite()
    {
        coordinates.clear();
        waterByHour.clear();
        initialWater.clear();
    }
};

std::vector<double> HAMSTAD_Benchmark2::coordinates{};
std::map<int, std::vector<double>> HAMSTAD_Benchmark2::waterByHour{};
std::vector<double> HAMSTAD_Benchmark2::initialWater{};

//! The initial state the report prescribes, w = 84,7687 kg/m3 at 95 %, is what the
//! engine's tabulated isotherm starts from, and all three profiles were captured.
TEST_F(HAMSTAD_Benchmark2, StartsFromThePrescribedState)
{
    ASSERT_EQ(3u, waterByHour.size());
    for(const double water : initialWater)
    {
        EXPECT_NEAR(Bench::initialWater, water, 1e-2);
    }
}

/////////////////////////////////////////////////////////////////////////////////////
/// MEASURED RESULTS, recorded 2026-09-15 (1 mm cells, ten-minute steps for ten hours
/// then hourly). Deviation from the analytical column, kg/m3, at every printed
/// position: the engine runs a shade dry behind the drying fronts and never ahead.
///
///   hour   worst |diff|   where            printed scatter of the seven codes
///    100      0,056       x = 0,04 m       up to 0,9 (CTH at x = 0,01)
///    300      0,024       x = 0,06-0,07 m  up to 0,3
///   1000      0,011       x = 0,08-0,10 m  up to 0,1
///
/// The faces sit on the printed values to 1e-4: the prescribed beta_p = 1e-3 s/m is,
/// as the analytical column itself shows, no surface resistance at all.
/////////////////////////////////////////////////////////////////////////////////////

//! The three printed profiles, every printed position. Tolerance roughly twice the
//! worst measured deviation, per the folder's convention.
TEST_F(HAMSTAD_Benchmark2, WaterContentProfilesMatchTheAnalyticalSolution)
{
    constexpr double tolerance{0.1};
    for(const auto & [position, hour100, hour300, hour1000] : analytical)
    {
        EXPECT_NEAR(hour100, valueAt(coordinates, waterByHour.at(100), position), tolerance)
          << "100 h, x = " << position << " m";
        EXPECT_NEAR(hour300, valueAt(coordinates, waterByHour.at(300), position), tolerance)
          << "300 h, x = " << position << " m";
        EXPECT_NEAR(hour1000, valueAt(coordinates, waterByHour.at(1000), position), tolerance)
          << "1000 h, x = " << position << " m";
    }
}
