#include <algorithm>
#include <fstream>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

#include "BeamBuilder.hxx"
#include "TestHelpers.hxx"
#include "TestMaterials.hxx"

// ThermSample_* family: real-world THERM sample files reproduced as direct HygroThermFEM
// runs, to be exercised (gtest_filter=ThermSample_*) BEFORE any THERM build consumes a new
// HygroThermFEM. Source of truth for the geometry/BC values: the .thmz sample noted below.
//
// The THERM sample "Stucco Wall - Moisture" (Samples\Stucco Wall - Moisture.thmz) as a
// 1D beam: Stucco 25.4 mm | Laminated panel 19.05 mm | Fiberglass Batts 76.2 mm |
// Gypsum Board Interior 25.4 mm, exterior film 10 C / RH 0.10 / hc 7 on the stucco side,
// interior film 40 C / RH 0.30 / hc 10 on the gypsum side ("use first timestep" BCs),
// dt = 360 s. This is the real-world case reported to blow up when started from a
// near-saturated initial humidity (0.99); the mild sample initial condition (0.10) is
// kept as the companion case.
namespace
{
    struct WallRunResult
    {
        bool completed = false;
        std::string error;
        std::vector<std::vector<double>> temperatures;
        std::vector<std::vector<double>> waterContents;
        std::vector<std::vector<double>> humidities;
    };

    WallRunResult runStuccoWall(const double initialHumidity,
                                const unsigned numberOfSteps,
                                const std::string & csvPath = {},
                                const double dTime = 360.0)
    {
        WallRunResult result;
        HygroThermFEM::MultiDomain multiDomain;

        const HygroThermFEM::State initialState({.temperature = 21.0,
                                                 .humidity = initialHumidity,
                                                 .pressure = 101325.0,
                                                 .liquidPercent = 1.0});

        const auto & stucco = multiDomain.materials().createSolidMaterial(TestHelper::Stucco());
        const auto & panel =
          multiDomain.materials().createSolidMaterial(TestHelper::LaminatedPanel());
        const auto & fiberglass =
          multiDomain.materials().createSolidMaterial(TestHelper::FiberglassBatts());
        const auto & gypsum =
          multiDomain.materials().createSolidMaterial(TestHelper::GypsumBoardInterior());

        TestHelper::BeamBuilder builder(multiDomain);
        builder.xStart(0.0)
          .height(0.05)
          .numElementsY(1)
          .state(initialState)
          .addSegment({.material = stucco.name(), .numElementsX = 4, .width = 0.0254})
          .addSegment({.material = panel.name(), .numElementsX = 3, .width = 0.01905})
          .addSegment({.material = fiberglass.name(), .numElementsX = 6, .width = 0.0762})
          .addSegment({.material = gypsum.name(), .numElementsX = 4, .width = 0.0254})
          .build();

        const HygroThermFEM::FixedBCHCCoefficients exterior{10.0, 7.0, 0.1};
        const HygroThermFEM::FixedBCHCCoefficients interior{40.0, 10.0, 0.3};
        builder.applyBC_FixedHc(TestHelper::BeamBuilder::Edge::Left, exterior);
        builder.applyBC_FixedHc(TestHelper::BeamBuilder::Edge::Right, interior);

        try
        {
            auto temperatures =
              multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
            auto humidities = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);

            result.temperatures.push_back(temperatures);
            result.humidities.push_back(humidities);
            result.waterContents.push_back(
              multiDomain.nodes().properties(HygroThermFEM::Variable::water));

            for(unsigned step = 0; step < numberOfSteps; ++step)
            {
                const auto solution =
                  multiDomain.transient(temperatures, humidities, dTime, step);
                temperatures = solution.temperature;
                humidities = solution.humidity;
                result.temperatures.push_back(temperatures);
                result.humidities.push_back(humidities);
                result.waterContents.push_back(solution.waterContent);
            }
            result.completed = true;
        }
        catch(const std::exception & exc)
        {
            result.error = exc.what();
        }

        if(!csvPath.empty())
        {
            std::ofstream csv(csvPath);
            csv << "step,field";
            for(size_t idx = 0; idx < result.temperatures.front().size(); ++idx)
            {
                csv << ",n" << idx;
            }
            csv << "\n";
            for(size_t step = 0; step < result.temperatures.size(); ++step)
            {
                const auto writeRow = [&](const char * field, const std::vector<double> & row) {
                    csv << step << "," << field;
                    for(const auto val : row)
                    {
                        csv << "," << val;
                    }
                    csv << "\n";
                };
                writeRow("T", result.temperatures[step]);
                writeRow("phi", result.humidities[step]);
                writeRow("w", result.waterContents[step]);
            }
        }
        return result;
    }

    double maxAbsValue(const std::vector<std::vector<double>> & values)
    {
        double result = 0.0;
        for(const auto & row : values)
        {
            for(const auto val : row)
            {
                result = (std::max)(result, std::abs(val));
            }
        }
        return result;
    }
}   // namespace

TEST(ThermSample_StuccoWall, SampleInitialConditions)
{
    // The stock sample: initial 21 C / RH 0.10. Must run to completion with bounded,
    // physical fields.
    const auto result = runStuccoWall(0.1, 100);
    ASSERT_TRUE(result.completed) << result.error;

    for(const auto & temps : result.temperatures)
    {
        for(const auto val : temps)
        {
            EXPECT_GT(val, 10.0 - 2.0);
            EXPECT_LT(val, 40.0 + 2.0);
        }
    }
    for(const auto & phis : result.humidities)
    {
        for(const auto val : phis)
        {
            EXPECT_GE(val, 0.0);
            EXPECT_LE(val, 1.0);
        }
    }
}

TEST(ThermSample_StuccoWall, NearSaturatedNoLatentHeat)
{
    // Discriminator for the latent-heat sign defect: the same near-saturated start with
    // the heat-of-evaporation term excluded. If this stays bounded while the full-physics
    // variant self-heats above both boundary temperatures, the interior latent term is
    // the driver (heatOfEvaporation() dropped the minus sign of Theoretical Model eq. 59).
    HygroThermFEM::SimulationProperties::Instance().setCalculationParameters(
      false, true, false, false, false);
    const auto result = runStuccoWall(0.99, 100);
    HygroThermFEM::SimulationProperties::Instance().resetCalculationParameters();

    ASSERT_TRUE(result.completed) << result.error;
    std::cout << "[StuccoWall 0.99, no latent] max|T| over run = "
              << maxAbsValue(result.temperatures) << "\n";
    for(const auto & temps : result.temperatures)
    {
        for(const auto val : temps)
        {
            EXPECT_GT(val, 10.0 - 2.0);
            EXPECT_LT(val, 40.0 + 2.0);
        }
    }
}

TEST(ThermSample_StuccoWall, DISABLED_NearSaturatedRefinedDt)
{
    // dt-refinement probe (diagnostic, run explicitly with --gtest_also_run_disabled_tests):
    // same case at dt = 90 s. If the sub-freezing surface dip and the +latent overshoot
    // shrink toward the wet-bulb/condensation bounds, they are big-dt transients.
    const auto result = runStuccoWall(0.99, 400, "/tmp/htf_stucco_wall_099_dt90.csv", 90.0);
    ASSERT_TRUE(result.completed) << result.error;
    double minTemp = 1e9;
    double maxTemp = -1e9;
    for(const auto & temps : result.temperatures)
    {
        for(const auto val : temps)
        {
            minTemp = (std::min)(minTemp, val);
            maxTemp = (std::max)(maxTemp, val);
        }
    }
    std::cout << "[StuccoWall 0.99, dt=90] T range over run = [" << minTemp << ", "
              << maxTemp << "]" << std::endl;
}

TEST(ThermSample_StuccoWall, NearSaturatedInitialConditions)
{
    // The formerly reported blow-up: the wall started from RH 0.99. Before the latent-term
    // corrections (sign of h_lg per Theoretical Model eq. 59; flux-consistent vapour
    // potential c_sat*phi) this self-heated to the 1000-clamp; now temperatures stay
    // within the latent-effect envelope: evaporative cooling can reach the exterior air's
    // wet-bulb temperature (~2 C for 10 C / RH 0.10, briefly overshot at this coarse
    // dt = 360 s) and condensation warms a few degrees above the warm boundary. The
    // DISABLED_NearSaturatedRefinedDt probe documents convergence to [2.2, 41.3] C at
    // dt = 90 s. NOTE: transients crossing 0 C make the (currently disabled) freezing
    // model (D5) relevant for this wall in colder climates.
    const auto result = runStuccoWall(0.99, 100, "/tmp/htf_stucco_wall_099.csv");
    ASSERT_TRUE(result.completed) << result.error;

    std::cout << "[StuccoWall 0.99] max|T| over run = " << maxAbsValue(result.temperatures)
              << ", final T range = ["
              << *std::min_element(result.temperatures.back().begin(),
                                   result.temperatures.back().end())
              << ", "
              << *std::max_element(result.temperatures.back().begin(),
                                   result.temperatures.back().end())
              << "]\n";

    for(const auto & temps : result.temperatures)
    {
        for(const auto val : temps)
        {
            EXPECT_GT(val, -5.0);   // wet-bulb floor plus coarse-dt transient margin
            EXPECT_LT(val, 45.0);   // condensation-warming ceiling plus margin
        }
    }
    for(const auto & phis : result.humidities)
    {
        for(const auto val : phis)
        {
            EXPECT_GE(val, 0.0);
            EXPECT_LE(val, 1.0);
        }
    }
}
