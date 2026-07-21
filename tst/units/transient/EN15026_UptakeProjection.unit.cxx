#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "DumpCsv.hxx"
#include "EN15026Material.hxx"
#include "HygroThermFEM2D.hxx"
#include "SlabCreator.hxx"
#include "TestHelpers.hxx"

/////////////////////////////////////////////////////////////////////////////////////
/// EN 15026:2007 Annex A -- moisture uptake in a semi-infinite region, run as the
/// SCALAR-MU PROJECTION of the benchmark (see EN15026Material.hxx: the annex's
/// mu(w) spans 212..875 over the benchmark's own moisture range and the material
/// model carries one value, so the standard's Tables A.1/A.2 are NOT the pass
/// criterion here -- the validation book shows no scalar mu satisfies them). The
/// engine is instead compared against the 1D reference solver running the SAME
/// projection (hygrothermfem_python, case en15026_uptake_projection), which
/// isolates engine discretization and coupling from the representation question.
///
/// Setup per the annex: uniform initial state 20 C / phi = 0.5; at t = 0 the
/// x = 0 surface steps to 30 C / phi = 0.95 (both fields Dirichlet); the far end
/// is natural and far enough (30 m, geometric grading) to stay undisturbed over
/// the simulated week. Moisture-dependent lambda(w) is active via the tabular
/// 2D-conductivity branch. Bottom-row profiles at day 1 and day 7 are dumped
/// for the book's engine dataset.
///
/// This configuration -- both fields pinned on one edge of a coupled run --
/// exposed two penalty-boundary defects (fixed 2026-07-21, found by exactly this
/// test): TemperatureBC carried the physical vapour-flux energy term scaled by
/// its 1e18 numerical film coefficient (a phantom latent sink that dragged the
/// pinned temperature off its value), and MoistureBCFixedHumidity pinned the
/// vapour CONTENT c_sat(T_air) phi_air rather than the humidity (off by the
/// saturation ratio until the surface reaches the boundary temperature). Either
/// defect stalls the nonlinear iteration against the physical clamp; with both
/// fixes the week runs in seconds.
/////////////////////////////////////////////////////////////////////////////////////

namespace
{
    //! Geometric mesh from the wetted surface: dx_first at x = 0 growing by
    //! `growth` per element until `minLength` is covered -- the same rule as the
    //! reference solver's Grid1D.graded (141 elements here).
    std::vector<double> gradedCoordinates(double dxFirst, double growth, double minLength)
    {
        std::vector<double> coords{0.0};
        double width{dxFirst};
        while(coords.back() < minLength)
        {
            coords.push_back(coords.back() + width);
            width *= growth;
        }
        return coords;
    }
}   // namespace

TEST(EN15026_UptakeProjection, SevenDaysScalarMu)
{
    SCOPED_TRACE("Begin Test: EN 15026 uptake, scalar-mu projection, 7 days.");

    // Enable the moisture/temperature-dependent conductivity branch (lambda(w));
    // every physical term stays included.
    HygroThermFEM::SimulationProperties::Instance().setCalculationParameters(
      false, false, false, false, true);
    HygroThermFEM::SimulationProperties::Instance().setIterationParameters(1.0, 1e-9, 500);

    {
        HygroThermFEM::MultiDomain multiDomain(
          {.performThermal = true, .performMoisture = true});

        constexpr double initialTemperature = 20.0;
        constexpr double initialHumidity = 0.5;
        constexpr double surfaceTemperature = 30.0;
        constexpr double surfaceHumidity = 0.95;

        constexpr HygroThermFEM::State state({
            .temperature = initialTemperature,
            .humidity = initialHumidity,
            .pressure = 101325.0,
            .liquidPercent = 1.0
        });

        const auto & material =
          multiDomain.materials().createSolidMaterial(TestHelper::EN15026AnnexA());

        const auto coords = gradedCoordinates(5.0e-4, 1.06, 30.0);
        const auto nColumns = coords.size();

        TestHelper::SlabBuilder(multiDomain)
          .gridXCoordinates(coords)
          .height(0.05)
          .material(material.name())
          .state(state)
          .build();

        // Both fields pinned at the x = 0 edge (1-based node indices: bottom 1, top 2).
        multiDomain.thermal().createBC_FixedTemperature(1, 2, surfaceTemperature);
        multiDomain.moisture().createBC_FixedHumidity(
          1, 2, HygroThermFEM::TemperatureAndHumidity{surfaceTemperature, surfaceHumidity});

        // Backward-Euler schedule refined at the initial step; block boundaries
        // land exactly on day 1 (step 258) and day 7 (step 330).
        const std::vector<std::pair<double, int>> timeBlocks{
          {30.0, 120}, {600.0, 138}, {7200.0, 72}};

        auto temperatures =
          multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
        auto humidities = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);

        std::vector<std::vector<double>> humidityAtDays;
        std::vector<std::vector<double>> temperatureAtDays;
        int step{0};
        for(const auto & [dTime, nSteps] : timeBlocks)
        {
            for(int blockStep = 0; blockStep < nSteps; ++blockStep)
            {
                const auto solution =
                  multiDomain.transient(temperatures, humidities, dTime, step);
                temperatures = solution.temperature;
                humidities = solution.humidity;
                ++step;
                if(step == 258 || step == 330)
                {
                    humidityAtDays.push_back(humidities);
                    temperatureAtDays.push_back(temperatures);
                }
            }
        }
        ASSERT_EQ(2u, humidityAtDays.size());

        TestHelper::CsvDump humidityDump("en15026_uptake_projection_humidity.csv", nColumns);
        TestHelper::CsvDump temperatureDump("en15026_uptake_projection_temperature.csv",
                                            nColumns);
        for(std::size_t row = 0; row < humidityAtDays.size(); ++row)
        {
            humidityDump.addRow(row + 1,
                                TestHelper::bottomRow(humidityAtDays[row], nColumns, 2));
            temperatureDump.addRow(
              row + 1, TestHelper::bottomRow(temperatureAtDays[row], nColumns, 2));
        }

        // Sanity: the pins hold, the far end is undisturbed, and both day-7
        // profiles decay monotonically from the wetted surface.
        const auto humidityDay7 = TestHelper::bottomRow(humidityAtDays[1], nColumns, 2);
        const auto temperatureDay7 =
          TestHelper::bottomRow(temperatureAtDays[1], nColumns, 2);
        EXPECT_NEAR(surfaceHumidity, humidityDay7.front(), 1e-9);
        EXPECT_NEAR(surfaceTemperature, temperatureDay7.front(), 1e-9);
        EXPECT_NEAR(initialHumidity, humidityDay7.back(), 1e-9);
        EXPECT_NEAR(initialTemperature, temperatureDay7.back(), 1e-3);
        for(std::size_t col = 1; col < nColumns; ++col)
        {
            EXPECT_LE(humidityDay7[col], humidityDay7[col - 1] + 1e-12) << "column " << col;
            EXPECT_LE(temperatureDay7[col], temperatureDay7[col - 1] + 1e-9)
              << "column " << col;
        }
    }

    HygroThermFEM::SimulationProperties::Instance().resetIterationParameters();
    HygroThermFEM::SimulationProperties::Instance().resetCalculationParameters();
}
