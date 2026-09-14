#include <cmath>
#include <vector>

#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "SlabCreator.hxx"
#include "TestHelpers.hxx"
#include "TestMaterials.hxx"

/////////////////////////////////////////////////////////////////////////////////////
/// The vapour diffusion resistance factor in its two forms.
///
/// A material supplies the factor either as a single value or as a curve in water
/// content. Both are held by one class so the element assembly builds the same
/// concrete type either way; these tests pin the two properties that makes that safe,
/// without reference to any particular standard or material.
///
///   1. A FLAT curve must reproduce the scalar exactly. If it does not, the curve path
///      is not evaluating what it claims, and every result computed through it is
///      suspect. This is the equivalence that lets an existing material be re-expressed
///      as a curve with no change in answer.
///   2. A RISING curve must dry the material relative to a scalar pinned at the curve's
///      dry-range value. More resistance where the material is wet means less vapour
///      arrives, which is the entire point of carrying the curve.
///
/// EN15026_AnnexABenchmark is where the feature is validated against published values;
/// what is asserted here is the mechanism.
/////////////////////////////////////////////////////////////////////////////////////

namespace
{
    //! The linear, non-porous, vapour-only material: the resistance factor is the ONLY
    //! thing moving moisture, which is what these tests want to isolate. It is also what
    //! keeps them in milliseconds -- driving a real sorption curve up its saturation
    //! knee turns every step into a hard nonlinear solve and buys no extra confidence
    //! about whether a table is being read.
    HygroThermFEM::SolidMaterialParams vapourOnlyMaterial(const double resistanceFactor)
    {
        auto params = TestHelper::LinearSorption();
        params.diffusionResistanceFactor = resistanceFactor;
        return params;
    }

    //! Water content spans 0 to 1 kg/m3 for this material, since its sorption curve is
    //! the identity; resistance-factor curves below are tabulated over that range.
    constexpr double maxWaterContent{1.0};

    //! Marches a wetting slab and returns the bottom-row water content at the end.
    //! Everything except the resistance factor is held identical between runs.
    std::vector<double> finalWaterProfile(const HygroThermFEM::SolidMaterialParams & params)
    {
        HygroThermFEM::MultiDomain multiDomain(
          {.performThermal = false, .performMoisture = true});
        const auto & material = multiDomain.materials().createSolidMaterial(params);

        const std::vector<double> coords{
          0.0, 0.002, 0.004, 0.006, 0.008, 0.010, 0.014, 0.018, 0.024, 0.030, 0.040, 0.060};
        const auto nColumns = coords.size();

        TestHelper::SlabBuilder(multiDomain)
          .gridXCoordinates(coords)
          .height(0.02)
          .material(material.name())
          .state({.temperature = 20.0,
                  .humidity = 0.4,
                  .pressure = 101325.0,
                  .liquidPercent = 1.0})
          .build();

        multiDomain.moisture().createBC_FixedHumidity(
          1, 2, HygroThermFEM::TemperatureAndHumidity{20.0, 0.9});

        auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
        auto humidities = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);
        // Twelve hours is already enough for the front to reach the middle of the slab,
        // which is all these comparisons need.
        std::vector<double> water;
        for(int step = 0; step < 12; ++step)
        {
            const auto solution = multiDomain.transient(temperatures, humidities, 3600.0, step);
            temperatures = solution.temperature;
            humidities = solution.humidity;
            water = TestHelper::bottomRow(solution.waterContent, nColumns, 2);
        }
        return water;
    }
}   // namespace

TEST(MoistureDependentResistanceFactor, FlatCurveReproducesTheScalar)
{
    SCOPED_TRACE("Begin Test: a constant mu(w) curve must equal the scalar it repeats.");

    constexpr double resistanceFactor{12.0};

    const auto scalarParams = vapourOnlyMaterial(resistanceFactor);

    auto curveParams = scalarParams;
    // Same number at both ends, so the table can only ever return that number.
    curveParams.diffusionResistanceFactorMoistureDependent = {
      {0.0, resistanceFactor}, {maxWaterContent, resistanceFactor}};

    const auto scalarProfile = finalWaterProfile(scalarParams);
    const auto curveProfile = finalWaterProfile(curveParams);

    ASSERT_EQ(scalarProfile.size(), curveProfile.size());
    for(std::size_t column = 0u; column < scalarProfile.size(); ++column)
    {
        EXPECT_NEAR(scalarProfile[column], curveProfile[column], 1e-9)
          << "a flat curve disagrees with the scalar at column " << column;
    }
}

TEST(MoistureDependentResistanceFactor, RisingCurveResistsMoreThanItsDryValue)
{
    SCOPED_TRACE("Begin Test: mu rising with water content must slow the wetting.");

    constexpr double dryResistanceFactor{12.0};

    const auto scalarParams = vapourOnlyMaterial(dryResistanceFactor);

    auto curveParams = scalarParams;
    // Equal to the scalar when dry, four times it when wet.
    curveParams.diffusionResistanceFactorMoistureDependent = {
      {0.0, dryResistanceFactor}, {maxWaterContent, 4.0 * dryResistanceFactor}};

    const auto scalarProfile = finalWaterProfile(scalarParams);
    const auto curveProfile = finalWaterProfile(curveParams);

    ASSERT_EQ(scalarProfile.size(), curveProfile.size());

    // Behind the wetted face the curve must hold moisture back. The face itself is
    // pinned by the boundary condition and the far end has not been reached, so the
    // comparison is made where the front actually is.
    bool anyDrier{false};
    for(std::size_t column = 1u; column + 1u < scalarProfile.size(); ++column)
    {
        EXPECT_LE(curveProfile[column], scalarProfile[column] + 1e-9)
          << "the rising curve let MORE moisture in than its own dry value, at column "
          << column;
        anyDrier = anyDrier || curveProfile[column] < scalarProfile[column] - 1e-6;
    }
    EXPECT_TRUE(anyDrier) << "the curve changed nothing: it is probably not being read at all";
}

TEST(MoistureDependentResistanceFactor, CurveAloneSatisfiesTheMaterialDataCheck)
{
    SCOPED_TRACE("Begin Test: a material carrying only the curve is fully specified.");

    auto params = vapourOnlyMaterial(0.0);
    params.diffusionResistanceFactorMoistureDependent = {{0.0, 10.0}, {maxWaterContent, 40.0}};

    HygroThermFEM::MultiDomain multiDomain({.performThermal = false, .performMoisture = true});
    const auto & material = multiDomain.materials().createSolidMaterial(params);

    EXPECT_TRUE(material.hasDiffusionResistanceFactorMoistureDependent());
    EXPECT_TRUE(material.hasAnyDiffusionResistanceFactor());
}
