#include <gtest/gtest.h>

#include "EN15026Material.hxx"
#include "HAMSTADBenchmarks.hxx"

/////////////////////////////////////////////////////////////////////////////////////
/// The HAMSTAD material transcriptions against the values the package itself prints,
/// and against the one material two documents print: benchmark 1's load bearing
/// material is the material EN 15026:2007 Annex A adopted, and EN15026Material.hxx
/// transcribed it from the standard independently of HAMSTADBenchmarks.hxx. Where the
/// two helpers disagree, one of the transcriptions is wrong. Runs in microseconds; see
/// the folder README on why property checks live next to the solves they protect.
/////////////////////////////////////////////////////////////////////////////////////

namespace
{
    namespace Ham = TestHelper::HAMSTAD;
    namespace Annex = TestHelper::EN15026;
}   // namespace

//! Benchmark 2, report p. 20-21: the isotherm's printed anchors. 84,7687 kg/m3 is the
//! prescribed initial state at 95 %; 19,5356 and 30,5059 are the analytical column's
//! face values, the equilibria with 45 % and 65 % air.
TEST(HAMSTAD_Properties, Benchmark2IsothermReproducesThePrintedAnchors)
{
    EXPECT_NEAR(84.7687, Ham::Benchmark2::waterContent(0.95), 5e-4);
    EXPECT_NEAR(19.5356, Ham::Benchmark2::waterContent(0.45), 5e-4);
    EXPECT_NEAR(30.5059, Ham::Benchmark2::waterContent(0.65), 5e-4);
    for(const double humidity : {0.45, 0.65, 0.80, 0.95})
    {
        EXPECT_NEAR(
          humidity,
          Ham::Benchmark2::humidityFromWaterContent(Ham::Benchmark2::waterContent(humidity)),
          1e-12);
    }
}

//! Benchmark 1's load bearing material against EN15026Material.hxx, function by
//! function, across the range both benchmarks visit.
TEST(HAMSTAD_Properties, Benchmark1LoadBearingMaterialIsTheEN15026Material)
{
    const Ham::Material material{Ham::Benchmark1::loadBearing()};
    for(const double humidity : {0.3, 0.5, 0.7, 0.9, 0.95, 0.99, 0.999})
    {
        const double suction{Ham::suctionFromHumidity(humidity)};
        const double water{Ham::waterContent(material, suction)};
        EXPECT_NEAR(Annex::waterContent(humidity), water, 1e-9) << "phi = " << humidity;
        EXPECT_NEAR(Annex::vapourResistanceFactor(water)
                      / Ham::Benchmark1::loadBearing().dryResistanceFactor,
                    Ham::vapourDiffusivityInAir / material.dryResistanceFactor
                      / Ham::vapourDiffusivity(material, water),
                    1e-9)
          << "phi = " << humidity;
        EXPECT_NEAR(Annex::liquidConductivity(water),
                    material.liquidConductivity(water),
                    1e-9 * Annex::liquidConductivity(water))
          << "phi = " << humidity;
        EXPECT_NEAR(Annex::moistureDiffusivity(water),
                    Ham::moistureDiffusivity(material, suction),
                    1e-6 * Annex::moistureDiffusivity(water))
          << "phi = " << humidity;
        EXPECT_NEAR(
          Annex::thermalConductivity(water), Ham::thermalConductivity(material, water), 1e-12)
          << "phi = " << humidity;
    }
    EXPECT_NEAR(Annex::dryVolumetricHeatCapacity, material.dryVolumetricHeatCapacity, 1e-6);
    const auto [density, heatCapacity] = Ham::capacitySplit(material.dryVolumetricHeatCapacity);
    EXPECT_NEAR(Annex::densitySplit, density, 1.0);
}

//! Benchmark 1, report p. 11: the printed inverse retention curves return the
//! prescribed initial water contents through the forward curves.
TEST(HAMSTAD_Properties, Benchmark1InverseRetentionCurvesRoundTrip)
{
    const Ham::Material loadBearing{Ham::Benchmark1::loadBearing()};
    const Ham::Material insulation{Ham::Benchmark1::insulation()};
    EXPECT_NEAR(
      145.0, Ham::waterContent(loadBearing, Ham::Benchmark1::loadBearingSuction(145.0)), 1e-9);
    EXPECT_NEAR(
      0.065, Ham::waterContent(insulation, Ham::Benchmark1::insulationSuction(0.065)), 1e-12);
}

//! Benchmark 5, report p. 47: at zero suction every retention curve returns w_sat
//! (the mode weights sum to one), and the analytic slope matches a central difference.
TEST(HAMSTAD_Properties, Benchmark5RetentionCurvesAreConsistent)
{
    for(const auto & material :
        {Ham::Benchmark5::brick(), Ham::Benchmark5::mortar(), Ham::Benchmark5::insulation()})
    {
        EXPECT_NEAR(material.freeSaturation, Ham::waterContent(material, 0.0), 1e-9)
          << material.name;
        for(const double suction : {1e3, 1e5, 1e7})
        {
            const double step{suction * 1e-4};
            const double central{(Ham::waterContent(material, suction + step)
                                  - Ham::waterContent(material, suction - step))
                                 / (2.0 * step)};
            const double analytic{
              Ham::retentionSlope(suction, material.freeSaturation, material.retention)};
            EXPECT_NEAR(central, analytic, 1e-5 * std::abs(central) + 1e-12)
              << material.name << " at P_suc = " << suction;
        }
    }
}
