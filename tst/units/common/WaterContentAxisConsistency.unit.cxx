#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

/////////////////////////////////////////////////////////////////////////////////////
/// The water-content axis check on a material's moisture tables.
///
/// Three tables are keyed by water content -- lambda(w), D_w(w) and mu(w) -- and each
/// is read at whatever water content the sorption curve returns, so its axis should
/// reach the curve's maximum. The check is advisory and lives outside the calculation
/// path: nothing in the engine calls it, a short table still evaluates through
/// extrapolation, and what these tests pin is that the report names the right table
/// with the right numbers, stays silent for constant tables, and works on a material
/// assembled through the setters in the order the mediator uses.
/////////////////////////////////////////////////////////////////////////////////////

namespace
{
    constexpr double maximumWaterContent{100.0};

    //! Every table closed exactly at the sorption curve's maximum.
    HygroThermFEM::SolidMaterialParams consistentMaterial()
    {
        return {.name = "porous",
                .thermalConductivityDry = 0.5,
                .density = 1000.0,
                .porosity = 0.2,
                .heatCapacity = 850.0,
                .diffusionResistanceFactor = 10.0,
                .thermalConductivityMoistureDependent = {{0.0, 0.5}, {maximumWaterContent, 0.9}},
                .moistureDependentMeasurementTemperature = 0.0,
                .thermalConductivityTemperatureDependent = {{0.0, 0.5}, {100.0, 0.5}},
                .temperatureDependentMeasurementHumidity = 0.0,
                .liquidTransportCurve = {{0.0, 0.0}, {maximumWaterContent, 1.0e-9}},
                .sorptionCurve = {{0.0, 0.0}, {0.5, 20.0}, {1.0, maximumWaterContent}}};
    }

    std::vector<std::string> tableNames(const std::vector<HygroThermFEM::WaterContentAxisMismatch> & mismatches)
    {
        std::vector<std::string> names;
        for(const auto & mismatch : mismatches)
        {
            names.push_back(mismatch.tableName);
        }
        return names;
    }
}   // namespace

TEST(WaterContentAxisConsistency, ConsistentTablesReportNothing)
{
    const HygroThermFEM::SolidMaterial material{consistentMaterial()};
    EXPECT_TRUE(HygroThermFEM::checkWaterContentAxes(material).empty());
}

TEST(WaterContentAxisConsistency, ShortLiquidTableIsReportedWithBothMaxima)
{
    auto params = consistentMaterial();
    params.liquidTransportCurve = {{0.0, 0.0}, {99.9, 1.0e-9}};
    const HygroThermFEM::SolidMaterial material{params};

    const auto mismatches = HygroThermFEM::checkWaterContentAxes(material);
    ASSERT_EQ(1u, mismatches.size());
    const auto & mismatch = mismatches.front();
    EXPECT_EQ("porous", mismatch.materialName);
    EXPECT_EQ("liquid transport coefficient", mismatch.tableName);
    EXPECT_DOUBLE_EQ(99.9, mismatch.tableMaximum);
    EXPECT_DOUBLE_EQ(maximumWaterContent, mismatch.sorptionMaximum);

    const auto message = mismatch.message();
    EXPECT_NE(std::string::npos, message.find("porous"));
    EXPECT_NE(std::string::npos, message.find("99.9"));
    EXPECT_NE(std::string::npos, message.find("100"));
}

TEST(WaterContentAxisConsistency, ShortResistanceFactorCurveIsReported)
{
    auto params = consistentMaterial();
    params.diffusionResistanceFactorMoistureDependent = {{0.0, 10.0}, {80.0, 50.0}};
    const HygroThermFEM::SolidMaterial material{params};

    const auto mismatches = HygroThermFEM::checkWaterContentAxes(material);
    ASSERT_EQ(1u, mismatches.size());
    EXPECT_EQ("moisture-dependent vapor resistance factor", mismatches.front().tableName);
    EXPECT_DOUBLE_EQ(80.0, mismatches.front().tableMaximum);
}

TEST(WaterContentAxisConsistency, EveryShortTableIsReportedInTableOrder)
{
    auto params = consistentMaterial();
    params.thermalConductivityMoistureDependent = {{0.0, 0.5}, {90.0, 0.9}};
    params.liquidTransportCurve = {{0.0, 0.0}, {95.0, 1.0e-9}};
    params.diffusionResistanceFactorMoistureDependent = {{0.0, 10.0}, {80.0, 50.0}};
    const HygroThermFEM::SolidMaterial material{params};

    const std::vector<std::string> expected{"moisture-dependent thermal conductivity",
                                            "liquid transport coefficient",
                                            "moisture-dependent vapor resistance factor"};
    EXPECT_EQ(expected, tableNames(HygroThermFEM::checkWaterContentAxes(material)));
}

TEST(WaterContentAxisConsistency, ConstantTablesAreNeverReported)
{
    // A single point is what the library stores for a constant property; the engine
    // duplicates it, leaving both points at one water content. Valid at any w.
    auto params = consistentMaterial();
    params.thermalConductivityMoistureDependent = {{0.0, 0.5}};
    params.liquidTransportCurve = {{0.0, 0.0}};
    params.diffusionResistanceFactorMoistureDependent = {{0.0, 10.0}};
    const HygroThermFEM::SolidMaterial material{params};

    EXPECT_TRUE(HygroThermFEM::checkWaterContentAxes(material).empty());
}

TEST(WaterContentAxisConsistency, SetterBuiltMaterialIsCheckedOnItsFinalState)
{
    // The mediator's order: curves first, sorption last. The check must see the
    // finished material, not whichever table happened to be set first.
    HygroThermFEM::Materials pool;
    auto & material = pool.createSolidMaterial("from setters");
    material.setLiquidTransportationCurve({{0.0, 0.0}, {99.0, 1.0e-9}});
    EXPECT_TRUE(HygroThermFEM::checkWaterContentAxes(material).empty())
      << "without a sorption curve there is nothing to compare against";

    material.setSorptionCurve({{0.0, 0.0}, {1.0, maximumWaterContent}});

    const auto mismatches = HygroThermFEM::checkWaterContentAxes(pool);
    ASSERT_EQ(1u, mismatches.size());
    EXPECT_EQ("from setters", mismatches.front().materialName);
    EXPECT_EQ("liquid transport coefficient", mismatches.front().tableName);
}
