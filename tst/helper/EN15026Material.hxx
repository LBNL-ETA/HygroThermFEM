#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

#include "HygroThermFEM2D.hxx"

//! EN 15026:2007 (E), Annex A (normative) -- "Benchmark example: moisture uptake
//! in a semi-infinite region".
//!
//! Every constant and every formula in this file is printed in the standard, and
//! each one carries the section it comes from: the general data of A.2
//! "a) General data" (page 18) and the property functions of A.2 "b) Material
//! data" (page 19). Nothing is taken from a secondary source, and nothing is a
//! fitted or stored number -- the tables the engine consumes are BUILT by
//! evaluating these functions, so the file can be read against the document.
//!
//! EN15026_AnnexAProperties asserts this transcription against the standard's own
//! published values. The two places where the engine's material model cannot hold what
//! the annex specifies -- the scalar diffusion resistance factor, and the
//! humidity-keyed conductivity table -- are called out at the functions concerned.
namespace TestHelper::EN15026
{
    /////////////////////////////////////////////////////////////////////////////
    /// A.2 a) General data (page 18)
    /////////////////////////////////////////////////////////////////////////////

    inline constexpr double molarGasConstant{8.314};        //!< R_H2O as printed, J/(mol K)
    inline constexpr double molarMassWater{0.018};          //!< M_w, kg/mol
    inline constexpr double densityWater{1000.0};           //!< rho_w, kg/m3
    inline constexpr double referenceTemperature{293.15};   //!< T = T_ref, K

    //! The gas constant as the SUCTION relation needs it, R / M_w = 461.9 J/(kg K).
    //!
    //! The annex prints R_H2O in J/(mol K) but then writes R_H2O T rho_w with rho_w
    //! given in kg/m3, and that product is only a pressure if the gas constant is the
    //! SPECIFIC one -- equivalently, if rho_w is read as a molar density. Both
    //! readings give the same number. The vapour-diffusion formula further down uses
    //! the MOLAR constant with M_w alongside it, which is where the printed 8,314
    //! belongs; the annex uses one symbol for both. Reading it as molar in the
    //! suction relation puts the benchmark's initial state at w = 143,8 kg/m3 instead
    //! of 42,9, which Table A.1's own undisturbed band of 40,8 to 45,1 rules out
    //! outright -- asserted in EN15026_AnnexAProperties rather than argued here.
    inline constexpr double specificGasConstantVapour{molarGasConstant / molarMassWater};

    /////////////////////////////////////////////////////////////////////////////
    /// A.2 b) Material data (page 19)
    /////////////////////////////////////////////////////////////////////////////

    //! w_max. "Porosity: equal maximum point of moisture storage function", so this
    //! one number is both the retention curve's numerator and the porosity source.
    inline constexpr double freeSaturation{146.0};   //!< kg/m3

    inline constexpr double retentionFactor{8.0e-8};   //!< the printed 8 x 10^-8, 1/Pa
    inline constexpr double retentionExponent{1.6};
    inline constexpr double retentionPower{0.375};

    inline constexpr double vapourDiffusivityInAir{26.1e-6};   //!< the printed 26,1 x 10^-6, m2/s
    inline constexpr double dryResistanceFactor{200.0};        //!< the printed 200

    inline constexpr double dryThermalConductivity{1.5};             //!< W/(m K)
    inline constexpr double conductivityMoistureSlope{15.8e-3};      //!< the printed 15,8 / 1000
    inline constexpr double dryVolumetricHeatCapacity{1.824e6};      //!< rho_0 c_0, J/(m3 K)

    //! Coefficients of the printed liquid-conductivity polynomial in (w - 73),
    //! lowest power first: K = exp(a0 + a1 (w-73) + ... + a5 (w-73)^5).
    inline constexpr std::array<double, 6> liquidConductivityCoefficients{
      -39.2619, 0.0704, -1.7420e-4, -2.7953e-6, -1.1566e-7, 2.5969e-9};

    /////////////////////////////////////////////////////////////////////////////
    /// Water retention and sorption
    /////////////////////////////////////////////////////////////////////////////

    //! Water retention curve, the annex's first formula:
    //!     w = 146 / (1 + (8e-8 p_suc)^1,6)^0,375
    [[nodiscard]] inline double waterContentFromSuction(const double suction)
    {
        const double scaled{std::pow(retentionFactor * suction, retentionExponent)};
        return freeSaturation / std::pow(1.0 + scaled, retentionPower);
    }

    //! The printed inverse p_suc(w), the annex's second formula:
    //!     p_suc = 0,125e8 ((146/w)^(1/0,375) - 1)^0,625
    //!
    //! Its constants are the exact inverse of the forward curve's -- 0,125e8 = 1/8e-8
    //! and 0,625 = 1/1,6 -- which is why they are written that way here rather than
    //! repeated as literals. RoundTripsThroughTheSuctionCurve pins the identity.
    [[nodiscard]] inline double suctionFromWaterContent(const double water)
    {
        const double ratio{std::pow(freeSaturation / water, 1.0 / retentionPower) - 1.0};
        return std::pow(ratio, 1.0 / retentionExponent) / retentionFactor;
    }

    //! Kelvin equation, the standard's eq. (6): p_suc = -rho_w R_H2O T ln(phi).
    //!
    //! Taken at T_ref because the standard's assumption list (4.1) makes the moisture
    //! storage function temperature independent, which is what the annex supplies
    //! T = T_ref for.
    [[nodiscard]] inline double suctionFromHumidity(const double humidity)
    {
        return -densityWater * specificGasConstantVapour * referenceTemperature
               * std::log(humidity);
    }

    //! Sorption isotherm w(phi), the annex's third formula -- which is the retention
    //! curve with Kelvin substituted, exactly as the annex prints it side by side.
    //! At phi = 1 the suction is zero and this returns the free saturation exactly.
    [[nodiscard]] inline double waterContent(const double humidity)
    {
        return humidity > 0.0 ? waterContentFromSuction(suctionFromHumidity(humidity)) : 0.0;
    }

    //! The printed inverse isotherm, the annex's fourth formula:
    //!     phi = exp(-p_suc(w) / (R_H2O T rho_w))
    [[nodiscard]] inline double humidityFromWaterContent(const double water)
    {
        return std::exp(-suctionFromWaterContent(water)
                        / (specificGasConstantVapour * referenceTemperature * densityWater));
    }

    //! dp_suc/dw [Pa m3/kg], the analytic derivative of the printed inverse. With
    //! u(w) = (146/w)^(1/0,375) - 1 the annex's inverse is p_suc = u^0,625 / 8e-8, so
    //!
    //!     du/dw     = -(1/0,375) (146/w)^(1/0,375) / w
    //!     dp_suc/dw = (0,625 / 8e-8) u^(-0,375) du/dw
    //!
    //! Negative throughout: suction falls as water content rises. Singular at
    //! w = w_max, where u reaches zero -- the benchmark never gets there, and the
    //! liquid table below stops short of it for that reason.
    [[nodiscard]] inline double suctionDerivative(const double water)
    {
        const double inversePower{1.0 / retentionPower};
        const double ratio{std::pow(freeSaturation / water, inversePower)};
        const double slope{-inversePower * ratio / water};
        return std::pow(ratio - 1.0, 1.0 / retentionExponent - 1.0) * slope
               / (retentionFactor * retentionExponent);
    }

    /////////////////////////////////////////////////////////////////////////////
    /// Transport coefficients
    /////////////////////////////////////////////////////////////////////////////

    //! The saturation-deficit factor in the annex's vapour-diffusion formula:
    //!     f(w) = (1 - w/146) / (0,503 (1 - w/146)^2 + 0,497)
    [[nodiscard]] inline double saturationDeficitFactor(const double water)
    {
        const double deficit{1.0 - water / freeSaturation};
        return deficit / (0.503 * deficit * deficit + 0.497);
    }

    //! Diffusion resistance factor mu(w) = 200 / f(w) [-].
    //!
    //! The annex prints the vapour permeability with 26,1e-6 / 200 and the factor f(w)
    //! as separate pieces; reading it as "still air divided by mu" makes the printed
    //! 200 the DRY-RANGE value and f the moisture dependence. mu diverges at free
    //! saturation, where f reaches zero; over the benchmark's own range it runs from
    //! 212 at the initial state to 866 at the wetted surface. That span is the whole
    //! reason the benchmark needs a resistance factor curve rather than a single value.
    [[nodiscard]] inline double vapourResistanceFactor(const double water)
    {
        return dryResistanceFactor / saturationDeficitFactor(water);
    }

    //! Vapour permeability delta_p(w, T) [kg/(m s Pa)], the annex's formula in full:
    //!     delta_p = M_w/(R T) 26,1e-6/200 (1 - w/146) / (0,503 (1 - w/146)^2 + 0,497)
    //!
    //! R is the MOLAR gas constant here, paired with M_w in kg/mol -- that pairing is
    //! what makes the result kg/(m s Pa). The engine takes mu rather than delta_p, so
    //! this function exists to state the printed formula and to let the test check
    //! that mu really is the ratio the annex implies.
    [[nodiscard]] inline double vapourPermeability(const double water, const double temperature)
    {
        return molarMassWater / (molarGasConstant * temperature) * vapourDiffusivityInAir
               / vapourResistanceFactor(water);
    }

    //! Liquid water permeability K(w) [s/m], the annex's printed polynomial:
    //!     K = exp(-39,2619 + 0,0704 (w-73) - 1,7420e-4 (w-73)^2 - 2,7953e-6 (w-73)^3
    //!             - 1,1566e-7 (w-73)^4 + 2,5969e-9 (w-73)^5)
    [[nodiscard]] inline double liquidConductivity(const double water)
    {
        const double shifted{water - 73.0};
        double exponent{0.0};
        double power{1.0};
        for(const double coefficient : liquidConductivityCoefficients)
        {
            exponent += coefficient * power;
            power *= shifted;
        }
        return std::exp(exponent);
    }

    //! Moisture diffusivity D_w(w) = -K(w) dp_suc/dw [m2/s], the standard's eq. (19)
    //! rewrite of the liquid flux g_w = -K grad(p_suc) into the water-content
    //! gradient. Positive, because dp_suc/dw is negative. This is the quantity the
    //! engine's liquid transport curve stores, and the conversion
    //! EN15026_LiquidUptakeSimilarity checks on its own.
    [[nodiscard]] inline double moistureDiffusivity(const double water)
    {
        return -liquidConductivity(water) * suctionDerivative(water);
    }

    //! Thermal conductivity lambda(w) = 1,5 + (15,8 / 1000) w [W/(m K)].
    [[nodiscard]] inline double thermalConductivity(const double water)
    {
        return dryThermalConductivity + conductivityMoistureSlope * water;
    }

    /////////////////////////////////////////////////////////////////////////////
    /// Representation choices the annex does not itself prescribe
    /////////////////////////////////////////////////////////////////////////////

    //! The annex prints only the PRODUCT rho_0 c_0 = 1,824e6 J/(m3 K), while the
    //! engine's capacity form is (rho + w)(c + (w / rho_w) c_p,w) and so needs a
    //! split. That form adds an excess term on top of the standard's rho_0 c_0 + w
    //! c_w; minimising the excess subject to rho c = rho_0 c_0 gives
    //! rho = sqrt(rho_0 c_0 rho_w / c_p,w) = 660,3 kg/m3. The rounded 660 keeps the
    //! excess within 3,3 % of the total capacity at the benchmark's initial water
    //! content.
    inline constexpr double densitySplit{660.0};   //!< kg/m3

    //! J/(kg K)
    inline constexpr double heatCapacitySplit{dryVolumetricHeatCapacity / densitySplit};

    //! Humidity samples geometrically refined toward saturation: (1 - phi) runs
    //! geometrically from 0,999 down to 5e-4, so the spacing tightens exactly where
    //! the isotherm steepens. A TABULATION choice, not the standard's; how well it
    //! represents the closed form is measured in EN15026_AnnexAProperties.
    [[nodiscard]] inline std::vector<double> humidityGrid(const std::size_t points)
    {
        std::vector<double> grid;
        grid.reserve(points);
        for(std::size_t index = 0u; index < points; ++index)
        {
            constexpr double firstDeficit{0.999};
            constexpr double lastDeficit{5.0e-4};
            const double fraction{static_cast<double>(index)
                                  / static_cast<double>(points - 1u)};
            grid.push_back(1.0 - firstDeficit * std::pow(lastDeficit / firstDeficit, fraction));
        }
        return grid;
    }

    /////////////////////////////////////////////////////////////////////////////
    /// Tabulation into the engine's material
    /////////////////////////////////////////////////////////////////////////////

    //! Sorption isotherm as (phi, w), closed at both ends by the values the formula
    //! itself takes there: dry at phi = 0, free saturation at phi = 1.
    [[nodiscard]] inline std::vector<FenestrationCommon::point>
      sorptionTable(const std::vector<double> & humidities)
    {
        std::vector<FenestrationCommon::point> table{{0.0, 0.0}};
        table.reserve(humidities.size() + 2u);
        for(const double humidity : humidities)
        {
            table.push_back({humidity, waterContent(humidity)});
        }
        table.push_back({1.0, freeSaturation});
        return table;
    }

    //! lambda keyed by RELATIVE HUMIDITY rather than by water content, because the
    //! engine builds its moisture-dependent conductivity with Variable::humidity. The
    //! annex's linear lambda(w) therefore enters as lambda(w(phi)) on the sorption
    //! grid -- the same curve, sampled through the isotherm.
    [[nodiscard]] inline std::vector<FenestrationCommon::point>
      conductivityTable(const std::vector<double> & humidities)
    {
        std::vector<FenestrationCommon::point> table{{0.0, thermalConductivity(0.0)}};
        table.reserve(humidities.size() + 2u);
        for(const double humidity : humidities)
        {
            table.push_back({humidity, thermalConductivity(waterContent(humidity))});
        }
        table.push_back({1.0, thermalConductivity(freeSaturation)});
        return table;
    }

    //! The annex's resistance factor as (w, mu). This is the table the engine reads when
    //! the material is built with the annex's real moisture dependence rather than the
    //! dry-range scalar. No saturation endpoint: f(w) reaches zero at w_max so mu diverges
    //! there, and the benchmark never gets closer than the last grid humidity.
    [[nodiscard]] inline std::vector<FenestrationCommon::point>
      resistanceFactorTable(const std::vector<double> & humidities)
    {
        std::vector<FenestrationCommon::point> table;
        table.reserve(humidities.size());
        for(const double humidity : humidities)
        {
            const double water{waterContent(humidity)};
            table.push_back({water, vapourResistanceFactor(water)});
        }
        return table;
    }

    //! Liquid transport as (w, D_w). No saturation endpoint: dp_suc/dw is singular at
    //! w_max, so the curve stops at the last grid humidity.
    [[nodiscard]] inline std::vector<FenestrationCommon::point>
      liquidTransportTable(const std::vector<double> & humidities)
    {
        std::vector<FenestrationCommon::point> table;
        table.reserve(humidities.size());
        for(const double humidity : humidities)
        {
            const double water{waterContent(humidity)};
            table.push_back({water, moistureDiffusivity(water)});
        }
        return table;
    }

    //! The Annex A material as the engine consumes it -- the closed-form functions
    //! above, tabulated on one shared humidity grid.
    //!
    //! TWO FORMS OF THE RESISTANCE FACTOR. The annex's mu is moisture dependent,
    //! mu(w) = 200/f(w), running 212 to 866 over the benchmark's own range. The engine
    //! can now carry it either way, so this builder offers both: the annex's real curve,
    //! which is what EN15026_AnnexABenchmark needs to meet Table A.1, and the dry-range
    //! scalar, which is the PROJECTION the engine was limited to before and which
    //! EN15026_UptakeProjection still pins.
    //!
    //! \param scalarMu The one diffusion resistance factor this material carries when
    //!        \p moistureDependentMu is false. Defaults to the annex's printed dry-range
    //!        value. Passing something enormous suppresses vapour transport altogether,
    //!        reducing the model to pure liquid diffusion -- how
    //!        EN15026_LiquidUptakeSimilarity exercises the standard's K to D_w conversion
    //!        on its own.
    //! \param moistureDependentMu True builds the material the annex actually specifies,
    //!        carrying mu(w) as a table; the scalar is then ignored by the engine. This is
    //!        the form EN15026_AnnexABenchmark needs to satisfy Table A.1.
    //! \param tablePoints Humidity samples per property curve. 64 holds the isotherm
    //!        to 0,09 kg/m3 over the benchmark's range, measured by
    //!        EN15026_AnnexAProperties; the engine reads these tables at every Gauss
    //!        point, so the grid is kept short deliberately.
    [[nodiscard]] inline HygroThermFEM::SolidMaterialParams
      material(const double scalarMu = dryResistanceFactor,
               const bool moistureDependentMu = false,
               const std::size_t tablePoints = 64u)
    {
        const auto humidities{humidityGrid(tablePoints)};
        return {.name = moistureDependentMu ? "EN 15026 Annex A"
                                            : "EN 15026 Annex A (scalar-mu projection)",
                .thermalConductivityDry = dryThermalConductivity,
                .density = densitySplit,
                .porosity = freeSaturation / densityWater,
                .heatCapacity = heatCapacitySplit,
                .diffusionResistanceFactor = scalarMu,
                .diffusionResistanceFactorMoistureDependent =
                  moistureDependentMu ? resistanceFactorTable(humidities)
                                      : std::vector<FenestrationCommon::point>{},
                .thermalConductivityMoistureDependent = conductivityTable(humidities),
                .moistureDependentMeasurementTemperature = 0,
                .thermalConductivityTemperatureDependent = {{0.0, dryThermalConductivity},
                                                            {100.0, dryThermalConductivity}},
                .temperatureDependentMeasurementHumidity = 0,
                .liquidTransportCurve = liquidTransportTable(humidities),
                .sorptionCurve = sorptionTable(humidities)};
    }

    //! Geometric mesh from the wetted surface: dxFirst at x = 0 growing by `growth`
    //! per element until `minLength` is covered -- the same rule as the reference
    //! solver's Grid1D.graded, so the two sides discretise the benchmark identically.
    //! A solver-side choice; the annex prescribes no mesh.
    [[nodiscard]] inline std::vector<double>
      gradedCoordinates(const double dxFirst, const double growth, const double minLength)
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
}   // namespace TestHelper::EN15026
