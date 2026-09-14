#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include <gtest/gtest.h>

#include "EN15026Material.hxx"
#include "HygroThermFEM2D.hxx"

/////////////////////////////////////////////////////////////////////////////////////
/// EN 15026:2007 Annex A -- the benchmark MATERIAL, checked against the standard.
///
/// tst/helper/EN15026Material.hxx transcribes the annex's printed property functions
/// and builds the engine's tables by evaluating them. This file is what makes that
/// transcription auditable: every expectation below is a value the standard itself
/// publishes, or an identity the standard's own formulas must satisfy. Nothing here
/// runs a solve; EN15026_AnnexABenchmark does that, against Tables A.1 and A.2.
///
/// Why this is worth its own file. The annex prints four formulas for the storage
/// function (retention curve, its inverse, the isotherm, its inverse), each with
/// constants that repeat in different forms, and it prints R_H2O in units that do not
/// balance in the relation that uses it. A transcription can be wrong in ways that
/// still run and still look plausible, and the benchmark itself is a poor detector:
/// it takes a coupled week-long solve to reach a number, and the engine cannot pass
/// the standard's Tables A.1 and A.2 anyway while its diffusion resistance factor is
/// a single scalar. These checks run in microseconds and fail loudly.
/////////////////////////////////////////////////////////////////////////////////////

namespace
{
    namespace Annex = TestHelper::EN15026;

    //! The benchmark's own two states, from A.2 "Problem description": the material
    //! starts in equilibrium at phi = 0,5 and the surface steps to phi = 0,95.
    constexpr double initialHumidity{0.5};
    constexpr double surfaceHumidity{0.95};

    //! Table A.1, the entries where the benchmark is still UNDISTURBED, so the limits
    //! bracket the initial state rather than the transient: day 7 at x = 0,04 m and
    //! day 30 at x = 0,06 m both read 40,8 to 45,1 kg/m3. Any transcription whose
    //! isotherm puts w(0,5) outside this band cannot pass the benchmark at all.
    constexpr double undisturbedMinimum{40.8};
    constexpr double undisturbedMaximum{45.1};
}   // namespace

//! The annex prints the retention curve and its inverse separately, with constants
//! written in different forms (8e-8 against 0,125e8, and 1,6 against 0,625). If either
//! is mistranscribed the pair stops being a true inverse, and the same holds for the
//! isotherm and its printed inverse. Nothing else in the suite would notice.
TEST(EN15026_AnnexAProperties, PrintedInversesAreTrueInverses)
{
    SCOPED_TRACE("Begin Test: EN 15026 Annex A printed inverse relations.");

    double worstWater{0.0};
    double worstHumidity{0.0};
    for(int sample = 1; sample < 1000; ++sample)
    {
        const double humidity{static_cast<double>(sample) / 1000.0};
        const double water{Annex::waterContent(humidity)};
        worstWater = (std::max)(
          worstWater,
          std::abs(Annex::waterContentFromSuction(Annex::suctionFromWaterContent(water)) - water));
        worstHumidity =
          (std::max)(worstHumidity, std::abs(Annex::humidityFromWaterContent(water) - humidity));
    }
    EXPECT_LT(worstWater, 1e-12) << "w -> p_suc -> w drifts by " << worstWater << " kg/m3";
    EXPECT_LT(worstHumidity, 1e-14) << "phi -> w -> phi drifts by " << worstHumidity;
}

//! The isotherm at the benchmark's two states, against the values the annex itself
//! shows. Figure A.1 plots every moisture profile rising to about 130 kg/m3 at the
//! wetted surface and settling near 40 to 45 in the undisturbed material, and
//! Table A.1 puts numbers on the second one.
TEST(EN15026_AnnexAProperties, IsothermReproducesTheBenchmarkStates)
{
    SCOPED_TRACE("Begin Test: EN 15026 Annex A sorption isotherm at the benchmark states.");

    const double initialWater{Annex::waterContent(initialHumidity)};
    const double surfaceWater{Annex::waterContent(surfaceHumidity)};

    EXPECT_GT(initialWater, undisturbedMinimum)
      << "w(0,5) = " << initialWater << " kg/m3 falls below Table A.1's undisturbed limit";
    EXPECT_LT(initialWater, undisturbedMaximum)
      << "w(0,5) = " << initialWater << " kg/m3 rises above Table A.1's undisturbed limit";

    // Figure A.1's surface value, read off the plot to the nearest few kg/m3.
    EXPECT_NEAR(130.0, surfaceWater, 2.0)
      << "w(0,95) = " << surfaceWater << " kg/m3 against Figure A.1's wetted surface";

    // At zero suction the curve must return the free saturation exactly, which is also
    // the annex's porosity statement.
    EXPECT_DOUBLE_EQ(Annex::freeSaturation, Annex::waterContent(1.0));
    EXPECT_DOUBLE_EQ(0.0, Annex::waterContent(0.0));
}

//! The annex prints R_H2O = 8,314 J/(mol K) and then uses it as R_H2O T rho_w with
//! rho_w in kg/m3, which is only a pressure for the SPECIFIC gas constant. Reading the
//! printed value literally is the single most likely transcription error in the whole
//! annex, so it gets its own test: the literal reading is excluded by Table A.1.
TEST(EN15026_AnnexAProperties, SuctionRelationNeedsTheSpecificGasConstant)
{
    SCOPED_TRACE("Begin Test: EN 15026 Annex A gas constant in the Kelvin relation.");

    EXPECT_NEAR(461.9, Annex::specificGasConstantVapour, 0.05);

    const double molarReading{Annex::waterContentFromSuction(
      -Annex::densityWater * Annex::molarGasConstant * Annex::referenceTemperature
      * std::log(initialHumidity))};
    EXPECT_GT(molarReading, undisturbedMaximum)
      << "the molar reading gives w(0,5) = " << molarReading
      << " kg/m3, which Table A.1 would have to bracket for it to be the right one";
}

//! Equation (19) converts the annex's K(w) into the diffusivity the engine stores, and
//! it consumes dp_suc/dw. The derivative is differentiated by hand in the helper, so it
//! is checked against a central difference of the printed inverse it came from.
TEST(EN15026_AnnexAProperties, SuctionDerivativeMatchesThePrintedInverse)
{
    SCOPED_TRACE("Begin Test: EN 15026 Annex A analytic dp_suc/dw.");

    double worstRelative{0.0};
    double worstWater{0.0};
    for(int sample = 1; sample < 200; ++sample)
    {
        const double water{Annex::freeSaturation * static_cast<double>(sample) / 200.0};
        const double step{water * 1e-6};
        const double difference{(Annex::suctionFromWaterContent(water + step)
                                 - Annex::suctionFromWaterContent(water - step))
                                / (2.0 * step)};
        const double analytic{Annex::suctionDerivative(water)};
        EXPECT_LT(analytic, 0.0) << "suction must fall as water content rises, at w = " << water;
        const double relative{std::abs(difference - analytic) / std::abs(analytic)};
        if(relative > worstRelative)
        {
            worstRelative = relative;
            worstWater = water;
        }
    }
    // Central differencing at this step size is itself good to about 1e-9 relative.
    EXPECT_LT(worstRelative, 1e-7)
      << "analytic dp_suc/dw disagrees with the printed inverse by " << worstRelative
      << " relative, worst at w = " << worstWater << " kg/m3";
}

//! The property that decides whether this benchmark can be met at all. The annex's
//! diffusion resistance factor moves by a factor of four across its own moisture range,
//! so a material model holding a single number cannot express it -- which is what kept
//! the engine off Table A.1 until the curve landed. Both ends are asserted so the span
//! stays a measured fact rather than a remark in a comment.
TEST(EN15026_AnnexAProperties, ResistanceFactorIsMoistureDependent)
{
    SCOPED_TRACE("Begin Test: EN 15026 Annex A mu(w) across the benchmark range.");

    EXPECT_DOUBLE_EQ(Annex::dryResistanceFactor, Annex::vapourResistanceFactor(0.0));
    EXPECT_NEAR(211.8, Annex::vapourResistanceFactor(Annex::waterContent(initialHumidity)), 0.1);
    EXPECT_NEAR(866.4, Annex::vapourResistanceFactor(Annex::waterContent(surfaceHumidity)), 0.1);

    // Both forms of the material are buildable. The default keeps the dry-range scalar,
    // which is the projection EN15026_UptakeProjection pins; asking for the moisture
    // dependence produces the curve EN15026_AnnexABenchmark runs on.
    const auto projection = Annex::material();
    EXPECT_DOUBLE_EQ(Annex::dryResistanceFactor, projection.diffusionResistanceFactor);
    EXPECT_TRUE(projection.diffusionResistanceFactorMoistureDependent.empty());

    const auto annexMaterial = Annex::material(Annex::dryResistanceFactor, true);
    ASSERT_FALSE(annexMaterial.diffusionResistanceFactorMoistureDependent.empty());
    const auto & curve = annexMaterial.diffusionResistanceFactorMoistureDependent;
    EXPECT_NEAR(211.8, curve.front().y, 60.0)
      << "the curve should start near the dry end of the benchmark's range";
    EXPECT_GT(curve.back().y, 800.0) << "the curve should reach the annex's wet-end factor";

    // delta_p is still air's permeability divided by mu, which is the reading of the
    // annex's formula the helper is built on.
    constexpr double temperature{293.15};
    const double stillAir{Annex::molarMassWater / (Annex::molarGasConstant * temperature)
                          * Annex::vapourDiffusivityInAir};
    const double water{Annex::waterContent(initialHumidity)};
    EXPECT_DOUBLE_EQ(stillAir / Annex::vapourResistanceFactor(water),
                     Annex::vapourPermeability(water, temperature));
}

//! The printed polynomial is written in (w - 73), so at w = 73 every term but the
//! constant vanishes and K must be exp(-39,2619) exactly. That is the one point where
//! the standard's own text pins the value without any arithmetic, and it catches a
//! transposed coefficient or a wrong shift.
TEST(EN15026_AnnexAProperties, LiquidConductivityMatchesThePrintedPolynomial)
{
    SCOPED_TRACE("Begin Test: EN 15026 Annex A liquid conductivity K(w).");

    EXPECT_DOUBLE_EQ(std::exp(-39.2619), Annex::liquidConductivity(73.0));

    // Equation (19)'s diffusivity must be positive and must rise with water content
    // across the benchmark's range; a sign slip in the derivative shows up here first.
    const double initialWater{Annex::waterContent(initialHumidity)};
    const double surfaceWater{Annex::waterContent(surfaceHumidity)};
    double previous{0.0};
    for(int sample = 0; sample <= 200; ++sample)
    {
        const double water{initialWater
                           + (surfaceWater - initialWater) * static_cast<double>(sample) / 200.0};
        const double diffusivity{Annex::moistureDiffusivity(water)};
        EXPECT_GT(diffusivity, previous) << "D_w is not rising at w = " << water << " kg/m3";
        previous = diffusivity;
    }
}

//! The two linear property statements, and the two representation choices the annex
//! leaves open, pinned where the material hands them to the engine.
TEST(EN15026_AnnexAProperties, LinearPropertiesAndTheCapacitySplit)
{
    SCOPED_TRACE("Begin Test: EN 15026 Annex A lambda(w), porosity and rho_0 c_0.");

    EXPECT_DOUBLE_EQ(1.5, Annex::thermalConductivity(0.0));
    EXPECT_DOUBLE_EQ(1.5 + 15.8 / 1000.0 * 146.0, Annex::thermalConductivity(146.0));

    const auto params = Annex::material();

    // "Porosity: equal maximum point of moisture storage function."
    EXPECT_DOUBLE_EQ(Annex::freeSaturation, params.porosity * Annex::densityWater);

    // The annex prints only the product; whatever split the helper chooses must keep it.
    EXPECT_DOUBLE_EQ(Annex::dryVolumetricHeatCapacity, params.density * params.heatCapacity);
    EXPECT_DOUBLE_EQ(1.5, params.thermalConductivityDry);
}

//! The engine reads tables, not formulas, so the last thing to establish is that the
//! tabulation represents the closed form. Sampling BETWEEN grid nodes over the
//! benchmark's own humidity range measures the interpolation error the engine actually
//! sees, which is otherwise an assumption.
TEST(EN15026_AnnexAProperties, TabulatedIsothermTracksTheClosedForm)
{
    SCOPED_TRACE("Begin Test: EN 15026 Annex A isotherm tabulation error.");

    const auto params = Annex::material();
    const auto & table = params.sorptionCurve;
    ASSERT_GT(table.size(), 2u);

    double worstError{0.0};
    double worstHumidity{0.0};
    for(int sample = 0; sample <= 2000; ++sample)
    {
        const double humidity{initialHumidity
                              + (surfaceHumidity - initialHumidity)
                                  * static_cast<double>(sample) / 2000.0};
        std::size_t upper{1u};
        while(upper + 1u < table.size() && table[upper].x < humidity)
        {
            ++upper;
        }
        const double leftHumidity{table[upper - 1u].x};
        const double leftWater{table[upper - 1u].y};
        const double rightHumidity{table[upper].x};
        const double rightWater{table[upper].y};
        const double fraction{(humidity - leftHumidity) / (rightHumidity - leftHumidity)};
        const double interpolated{leftWater + fraction * (rightWater - leftWater)};
        const double error{std::abs(interpolated - Annex::waterContent(humidity))};
        if(error > worstError)
        {
            worstError = error;
            worstHumidity = humidity;
        }
    }
    // Measured 0,089 kg/m3 at 64 points, against a profile spanning 43 to 129 kg/m3,
    // so the tabulation contributes about a tenth of a percent of the range. The
    // tolerance is set just above the measurement: a coarser grid or a changed grid
    // rule should have to be noticed here, not absorbed silently.
    EXPECT_LT(worstError, 0.12)
      << "linear interpolation of the sorption table misses the annex's isotherm by "
      << worstError << " kg/m3 at phi = " << worstHumidity;
}
