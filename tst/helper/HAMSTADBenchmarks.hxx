#pragma once

#include <cmath>
#include <cstddef>
#include <functional>
#include <string>
#include <utility>
#include <vector>

#include "HygroThermFEM2D.hxx"

//! HAMSTAD WP2 benchmark package -- the EU project "Determination of liquid water transfer
//! properties of porous building materials and development of numerical assessment
//! methods", GRD1-1999-20007. Two documents are transcribed here and every constant names
//! the page it is printed on:
//!
//!   [R]  Hagentoft, C-E. 2002. HAMSTAD - WP2 Modeling, Version 4. Report R-02:9,
//!        Department of Building Physics, Chalmers University of Technology. The
//!        mathematical model the benchmarks are defined against.
//!   [B]  Hagentoft, C-E. 2002. HAMSTAD - Final report: Methodology of HAM-modeling.
//!        Report R-02:8, Chalmers. Section 3 holds the five benchmark descriptions,
//!        each with its own page numbering carried below.
//!
//! Nothing here is taken from a secondary source and nothing is a fitted or stored
//! number: the tables the engine consumes are BUILT by evaluating the printed functions.
//! Where the engine's own representation differs from the report's (its temperature
//! dependent diffusion coefficient in air, its concentration-based surface coefficient,
//! its heat capacity form) the conversion is a named function with the reasoning at it.
namespace TestHelper::HAMSTAD
{
    /////////////////////////////////////////////////////////////////////////////
    /// General data, [B] p. 9 (benchmark 1) and p. 46 (benchmark 5)
    /////////////////////////////////////////////////////////////////////////////

    inline constexpr double gasConstant{8.314};             //!< R, J/(mol K)
    inline constexpr double molarMassWater{0.018};          //!< M_w, kg/mol
    inline constexpr double densityWater{1000.0};           //!< rho_w, kg/m3
    inline constexpr double referenceTemperature{293.15};   //!< T_ref, K
    inline constexpr double gravity{9.81};                  //!< g, m/s2

    //! R_v = R / M_w, the vapour gas constant the Kelvin relation uses, J/(kg K).
    inline constexpr double vapourGasConstant{gasConstant / molarMassWater};

    //! Diffusion coefficient of vapour in air the package fixes at T_ref, [B] p. 46:
    //! D_air,Schirmer = 2,31e-5 (T_ref / 273,16)^1,81 = 26,1e-6 m2/s. Every benchmark's
    //! vapour permeability is written with this one number.
    inline constexpr double vapourDiffusivityInAir{26.1e-6};

    /////////////////////////////////////////////////////////////////////////////
    /// Kelvin relation, [R] eq. (5)-(7) and [B] p. 46: P_suc = -rho_w R_v T ln(phi)
    /////////////////////////////////////////////////////////////////////////////

    //! Suction pressure [Pa] in equilibrium with a relative humidity, at T_ref. The
    //! package makes the sorption isotherm temperature independent ([R] section 3), which
    //! is what the reference temperature "for parameter transformations" exists for.
    [[nodiscard]] inline double suctionFromHumidity(const double humidity)
    {
        return -densityWater * vapourGasConstant * referenceTemperature * std::log(humidity);
    }

    //! The inverse: phi = exp(-P_suc / (rho_w R_v T_ref)).
    [[nodiscard]] inline double humidityFromSuction(const double suction)
    {
        return std::exp(-suction / (densityWater * vapourGasConstant * referenceTemperature));
    }

    /////////////////////////////////////////////////////////////////////////////
    /// Water retention, the multimodal van Genuchten form of [B] p. 46:
    ///     w(P_suc) = w_sat * sum_i k_i / (1 + (a_i P_suc)^n_i)^m_i,   m_i = 1 - 1/n_i
    /// Benchmark 1 prints its two materials in exactly this form with a single mode and
    /// the scale in 1/Pa ([B] p. 9); benchmark 5 prints the scale against the capillary
    /// suction HEIGHT h = P_suc / (rho_w g) in 1/m ([B] p. 46-47), converted at the
    /// material definition so that every mode below is keyed by pressure.
    /////////////////////////////////////////////////////////////////////////////

    struct RetentionMode
    {
        double weight;      //!< k_i [-]
        double scale;       //!< a_i [1/Pa]
        double exponentN;   //!< n_i [-]
        double exponentM;   //!< m_i [-]
    };

    //! w(P_suc) [kg/m3].
    [[nodiscard]] inline double retentionWater(const double suction,
                                               const double freeSaturation,
                                               const std::vector<RetentionMode> & modes)
    {
        double sum{0.0};
        for(const auto & [weight, scale, exponentN, exponentM] : modes)
        {
            sum += weight / std::pow(1.0 + std::pow(scale * suction, exponentN), exponentM);
        }
        return freeSaturation * sum;
    }

    //! dw/dP_suc [kg/(m3 Pa)], the analytic derivative of the form above. Negative
    //! throughout: water content falls as suction rises.
    [[nodiscard]] inline double retentionSlope(const double suction,
                                               const double freeSaturation,
                                               const std::vector<RetentionMode> & modes)
    {
        double sum{0.0};
        for(const auto & [weight, scale, exponentN, exponentM] : modes)
        {
            const double scaled{std::pow(scale * suction, exponentN)};
            const double outer{std::pow(1.0 + scaled, -exponentM - 1.0)};
            sum += -weight * exponentM * outer * exponentN * scaled / suction;
        }
        return freeSaturation * sum;
    }

    /////////////////////////////////////////////////////////////////////////////
    /// Vapour diffusion, [B] p. 9 and p. 46:
    ///     delta_p(w, T) = M_w / (R T) * 26,1e-6 / mu_dry * f(w)
    ///     f(w) = (1 - w/w_sat) / ((1 - p) (1 - w/w_sat)^2 + p)
    /// with p = 0,497 in benchmark 1 and p = 0,20 in benchmark 5.
    /////////////////////////////////////////////////////////////////////////////

    [[nodiscard]] inline double saturationDeficitFactor(const double water,
                                                        const double freeSaturation,
                                                        const double deficitFactorP)
    {
        const double deficit{1.0 - water / freeSaturation};
        return deficit / ((1.0 - deficitFactorP) * deficit * deficit + deficitFactorP);
    }

    /////////////////////////////////////////////////////////////////////////////
    /// One benchmark material as the report prints it
    /////////////////////////////////////////////////////////////////////////////

    //! Liquid water conductivity K(w) [s]; an empty function is a capillary non-active
    //! material (benchmark 1's insulation prints K = 0 s).
    using LiquidConductivity = std::function<double(double)>;

    //! K = exp(sum_i a_i (w - shift)^i), the polynomial of [B] p. 9 (shift 73 kg/m3).
    [[nodiscard]] inline LiquidConductivity
      exponentialPolynomialInWater(std::vector<double> coefficients, const double shift)
    {
        return [coefficients = std::move(coefficients), shift](const double water) {
            const double shifted{water - shift};
            double exponent{0.0};
            double power{1.0};
            for(const double coefficient : coefficients)
            {
                exponent += coefficient * power;
                power *= shifted;
            }
            return std::exp(exponent);
        };
    }

    //! K = exp(sum_i a_i (w / rho_w)^i), the polynomial of [B] p. 47.
    [[nodiscard]] inline LiquidConductivity
      exponentialPolynomialInVolumetricWater(std::vector<double> coefficients)
    {
        return [coefficients = std::move(coefficients)](const double water) {
            const double volumetric{water / densityWater};
            double exponent{0.0};
            double power{1.0};
            for(const double coefficient : coefficients)
            {
                exponent += coefficient * power;
                power *= volumetric;
            }
            return std::exp(exponent);
        };
    }

    struct Material
    {
        std::string name;
        double freeSaturation;   //!< w_sat [kg/m3]
        std::vector<RetentionMode> retention;
        double dryResistanceFactor;   //!< mu_dry [-]
        double deficitFactorP;        //!< p in f(w) [-]
        LiquidConductivity liquidConductivity;
        double dryThermalConductivity;      //!< lambda_dry [W/(m K)]
        double conductivityMoistureSlope;   //!< lambda = lambda_dry + slope * w, per kg/m3
        double dryVolumetricHeatCapacity;   //!< rho_0 c_0 [J/(m3 K)]
    };

    [[nodiscard]] inline double waterContent(const Material & material, const double suction)
    {
        return retentionWater(suction, material.freeSaturation, material.retention);
    }

    [[nodiscard]] inline double thermalConductivity(const Material & material, const double water)
    {
        return material.dryThermalConductivity + material.conductivityMoistureSlope * water;
    }

    //! D_vapour(w) = 26,1e-6 / mu_dry * f(w) [m2/s], the moisture dependent diffusion
    //! coefficient the report multiplies by M_w / (R T) to get delta_p ([B] p. 46).
    [[nodiscard]] inline double vapourDiffusivity(const Material & material, const double water)
    {
        return vapourDiffusivityInAir / material.dryResistanceFactor
               * saturationDeficitFactor(water, material.freeSaturation, material.deficitFactorP);
    }

    //! Moisture diffusivity D_w = -K dP_suc/dw [m2/s], the first expression of [R]
    //! eq. (37), which is what the engine's liquid transport curve stores. Written with
    //! the retention slope dw/dP_suc so that no curve inversion is needed: the tables
    //! below walk a suction grid and read w, phi and D_w off it side by side.
    [[nodiscard]] inline double moistureDiffusivity(const Material & material, const double suction)
    {
        const double water{waterContent(material, suction)};
        const double slope{retentionSlope(suction, material.freeSaturation, material.retention)};
        return -material.liquidConductivity(water) / slope;
    }

    /////////////////////////////////////////////////////////////////////////////
    /// Conversions into the engine's own conventions
    /////////////////////////////////////////////////////////////////////////////

    //! The engine's vapour flux is -(D_air(T) / mu) grad(c) on the vapour concentration
    //! c = p M_w / (R T), so its D_air / mu plays the role of the report's
    //! 26,1e-6 f(w) / mu_dry. The engine's D_air(T) is Hagentoft's (22,2 + 0,14 T) e-6,
    //! which reads 25,0e-6 at 20 C where the report pins 26,1e-6; the resistance factor
    //! handed to the engine is therefore D_air,engine(T_ref) / D_vapour(w), so that at
    //! T_ref the engine carries exactly the benchmark's permeability. Away from T_ref the
    //! engine keeps its temperature dependence where the report keeps none.
    [[nodiscard]] inline double engineResistanceFactor(const Material & material,
                                                       const double water)
    {
        const double engineDiffusivityInAir{
          HygroThermFEM::vaporDiffusionCoefficientAtTemperature(referenceTemperature - 273.15)};
        return engineDiffusivityInAir / vapourDiffusivity(material, water);
    }

    //! Surface vapour transfer. The report writes g = beta_p (p_air - p_surf) with beta_p
    //! in s/m ([R] eq. (45)); the engine writes g = beta_c (c_air - c_surf) with
    //! beta_c = h_c / (rho_air c_p,air), the Lewis relation, taking h_c from the
    //! FixedBCHCCoefficients it is given. With p = c R_v T the two agree at
    //! beta_c = beta_p R_v T_air, so the film coefficient that reproduces a prescribed
    //! beta_p is h = beta_p R_v T_air rho_air c_p,air. The product R_v T is read back from
    //! the engine's own saturation pair so the conversion cannot drift from it.
    [[nodiscard]] inline double filmCoefficientForVapourTransfer(const double betaP,
                                                                 const double airTemperature)
    {
        const double gasConstantTimesTemperature{
          HygroThermFEM::vaporPressureAtTemperature(airTemperature)
          / HygroThermFEM::saturationConcentrationAtTemperature(airTemperature)};
        return betaP * gasConstantTimesTemperature * Constants::Density_Air * Constants::Cp_Air;
    }

    //! The report prints the dry capacity as ONE product rho_0 c_0 while the engine's
    //! capacity is (rho + w)(c + (w / rho_w) c_p,w) and needs a split. That form carries an
    //! excess over the report's rho_0 c_0 + w c_w; minimising it subject to rho c = rho_0 c_0
    //! gives rho = sqrt(rho_0 c_0 rho_w / c_p,w), the choice EN15026Material.hxx made for
    //! the same material. Returned as (density, specific heat).
    [[nodiscard]] inline std::pair<double, double>
      capacitySplit(const double dryVolumetricHeatCapacity)
    {
        const double density{
          std::sqrt(dryVolumetricHeatCapacity * densityWater / Constants::Cp_Water)};
        return {density, dryVolumetricHeatCapacity / density};
    }

    /////////////////////////////////////////////////////////////////////////////
    /// Tabulation into the engine's material
    /////////////////////////////////////////////////////////////////////////////

    //! Suction samples, geometrically spaced from `minSuction` to `maxSuction` and returned
    //! DESCENDING, i.e. in order of rising humidity and water content -- the order every
    //! engine table wants. A tabulation choice, not the report's.
    [[nodiscard]] inline std::vector<double>
      suctionGrid(const std::size_t points, const double minSuction, const double maxSuction)
    {
        std::vector<double> grid;
        grid.reserve(points);
        for(std::size_t index = 0u; index < points; ++index)
        {
            const double fraction{static_cast<double>(index) / static_cast<double>(points - 1u)};
            grid.push_back(maxSuction * std::pow(minSuction / maxSuction, fraction));
        }
        return grid;
    }

    //! Sorption isotherm as (phi, w), closed at both ends by the values the retention
    //! curve itself takes there: dry at phi = 0, free saturation at phi = 1.
    [[nodiscard]] inline std::vector<FenestrationCommon::point>
      sorptionTable(const Material & material, const std::vector<double> & suctions)
    {
        std::vector<FenestrationCommon::point> table{{0.0, 0.0}};
        table.reserve(suctions.size() + 2u);
        for(const double suction : suctions)
        {
            table.emplace_back(humidityFromSuction(suction), waterContent(material, suction));
        }
        table.emplace_back(1.0, material.freeSaturation);
        return table;
    }

    //! lambda keyed by relative humidity, because the engine builds its moisture dependent
    //! conductivity with Variable::humidity: the report's linear lambda(w) enters as
    //! lambda(w(phi)) on the same suction grid.
    [[nodiscard]] inline std::vector<FenestrationCommon::point>
      conductivityTable(const Material & material, const std::vector<double> & suctions)
    {
        std::vector<FenestrationCommon::point> table{{0.0, material.dryThermalConductivity}};
        table.reserve(suctions.size() + 2u);
        for(const double suction : suctions)
        {
            table.emplace_back(humidityFromSuction(suction),
                               thermalConductivity(material, waterContent(material, suction)));
        }
        table.emplace_back(1.0, thermalConductivity(material, material.freeSaturation));
        return table;
    }

    //! The engine's resistance factor as (w, mu). No saturation endpoint: f(w) reaches
    //! zero at w_sat so mu diverges there.
    [[nodiscard]] inline std::vector<FenestrationCommon::point>
      resistanceFactorTable(const Material & material, const std::vector<double> & suctions)
    {
        std::vector<FenestrationCommon::point> table;
        table.reserve(suctions.size());
        for(const double suction : suctions)
        {
            const double water{waterContent(material, suction)};
            table.emplace_back(water, engineResistanceFactor(material, water));
        }
        return table;
    }

    //! Liquid transport as (w, D_w). A capillary non-active material (the report's
    //! K = 0 s) enters as a ZERO curve rather than no curve, the convention of the suite's
    //! other test materials: the engine's material constructor reads the liquid table
    //! unconditionally. No saturation endpoint otherwise: dw/dP_suc vanishes there and
    //! D_w diverges.
    [[nodiscard]] inline std::vector<FenestrationCommon::point>
      liquidTransportTable(const Material & material, const std::vector<double> & suctions)
    {
        std::vector<FenestrationCommon::point> table;
        if(!material.liquidConductivity)
        {
            table.emplace_back(0.0, 0.0);
            table.emplace_back(material.freeSaturation, 0.0);
            return table;
        }
        table.reserve(suctions.size());
        for(const double suction : suctions)
        {
            table.emplace_back(waterContent(material, suction),
                               moistureDiffusivity(material, suction));
        }
        return table;
    }

    //! The material as the engine consumes it, every property tabulated on one shared
    //! suction grid. 1e1 Pa is phi = 0,99999993 and 1e9 Pa is phi = 6e-4, which brackets
    //! every state the benchmarks visit including interstitial condensation.
    [[nodiscard]] inline HygroThermFEM::SolidMaterialParams
      solidMaterial(const Material & material, const std::size_t tablePoints = 100u)
    {
        const auto suctions{suctionGrid(tablePoints, 1.0e1, 1.0e9)};
        const auto [density, heatCapacity] = capacitySplit(material.dryVolumetricHeatCapacity);
        return {
          .name = material.name,
          .thermalConductivityDry = material.dryThermalConductivity,
          .density = density,
          .porosity = material.freeSaturation / densityWater,
          .heatCapacity = heatCapacity,
          .diffusionResistanceFactor = material.dryResistanceFactor,
          .diffusionResistanceFactorMoistureDependent = resistanceFactorTable(material, suctions),
          .thermalConductivityMoistureDependent = conductivityTable(material, suctions),
          .moistureDependentMeasurementTemperature = 0,
          .thermalConductivityTemperatureDependent = {{0.0, material.dryThermalConductivity},
                                                      {100.0, material.dryThermalConductivity}},
          .temperatureDependentMeasurementHumidity = 0,
          .liquidTransportCurve = liquidTransportTable(material, suctions),
          .sorptionCurve = sorptionTable(material, suctions)};
    }

    /////////////////////////////////////////////////////////////////////////////
    /// Benchmark 1, "Insulated roof", [B] p. 8-12
    /////////////////////////////////////////////////////////////////////////////
    namespace Benchmark1
    {
        //! A. Load bearing material, [B] p. 9. This is the material EN 15026:2007 Annex A
        //! later adopted verbatim; EN15026Material.hxx transcribes the same functions from
        //! the standard, and HAMSTAD_Benchmark1Properties checks the two agree.
        [[nodiscard]] inline Material loadBearing()
        {
            return {
              .name = "HAMSTAD BM1 load bearing",
              .freeSaturation = 146.0,
              .retention = {{.weight = 1.0, .scale = 8.0e-8, .exponentN = 1.6, .exponentM = 0.375}},
              .dryResistanceFactor = 200.0,
              .deficitFactorP = 0.497,
              .liquidConductivity = exponentialPolynomialInWater(
                {-39.2619, 0.0704, -1.7420e-4, -2.7953e-6, -1.1566e-7, 2.5969e-9}, 73.0),
              .dryThermalConductivity = 1.5,
              .conductivityMoistureSlope = 15.8 / 1000.0,
              .dryVolumetricHeatCapacity = 1.824e6};
        }

        //! B. Insulation material, [B] p. 9-10. Capillary non-active: K = 0 s.
        [[nodiscard]] inline Material insulation()
        {
            return {
              .name = "HAMSTAD BM1 insulation",
              .freeSaturation = 900.0,
              .retention = {{.weight = 1.0, .scale = 2.0e-4, .exponentN = 2.0, .exponentM = 0.5}},
              .dryResistanceFactor = 9.6,
              .deficitFactorP = 0.497,
              .liquidConductivity = {},
              .dryThermalConductivity = 0.033,
              .conductivityMoistureSlope = 0.59 / 1000.0,
              .dryVolumetricHeatCapacity = 0.0739e6};
        }

        //! The printed inverse retention curves, [B] p. 9: P_suc(w) for each material, used
        //! to turn the prescribed initial water contents into the humidity the engine solves
        //! for. The constants are the exact inverses of the forward curves' and are written
        //! that way.
        [[nodiscard]] inline double loadBearingSuction(const double water)
        {
            const double ratio{std::pow(146.0 / water, 1.0 / 0.375) - 1.0};
            return std::pow(ratio, 0.625) / 8.0e-8;
        }

        [[nodiscard]] inline double insulationSuction(const double water)
        {
            const double ratio{std::pow(900.0 / water, 2.0) - 1.0};
            return std::pow(ratio, 0.5) / 2.0e-4;
        }

        inline constexpr double loadBearingThickness{0.1};        //!< m, [B] p. 10
        inline constexpr double insulationThickness{0.05};        //!< m
        inline constexpr double exteriorFilmCoefficient{25.0};    //!< alpha_e,e W/(m2 K), [B] p. 11
        inline constexpr double interiorFilmCoefficient{7.0};     //!< alpha_e,i
        inline constexpr double exteriorVapourCoefficient{0.0};   //!< beta_p,e s/m: sealed
        inline constexpr double interiorVapourCoefficient{2.0e-8};   //!< beta_p,i s/m
        inline constexpr double initialLoadBearingWater{145.0};      //!< kg/m3, [B] p. 11
        inline constexpr double initialInsulationWater{0.065};       //!< kg/m3
        inline constexpr double initialTemperature{10.0};            //!< C
    }   // namespace Benchmark1

    /////////////////////////////////////////////////////////////////////////////
    /// Benchmark 2, "Analytical benchmark", [B] p. 19-21
    /////////////////////////////////////////////////////////////////////////////
    namespace Benchmark2
    {
        //! Sorption isotherm of material A, [B] p. 20:
        //!     w = 116 / (1 - ln(phi) / 0,118)^0,869
        [[nodiscard]] inline double waterContent(const double humidity)
        {
            return 116.0 / std::pow(1.0 - std::log(humidity) / 0.118, 0.869);
        }

        //! Its printed inverse: phi = exp(0,118 (1 - (116 / w)^(1/0,869))).
        [[nodiscard]] inline double humidityFromWaterContent(const double water)
        {
            return std::exp(0.118 * (1.0 - std::pow(116.0 / water, 1.0 / 0.869)));
        }

        inline constexpr double vapourPermeability{1.0e-15};    //!< delta_p, s (constant)
        inline constexpr double moistureDiffusivity{6.0e-10};   //!< D_w, m2/s (constant)
        inline constexpr double thermalConductivity{0.15};      //!< W/(m K)
        inline constexpr double heatCapacity{800.0};            //!< c_p, J/(kg K)
        inline constexpr double density{525.0};                 //!< rho_0, kg/m3
        inline constexpr double thickness{0.2};                 //!< m
        inline constexpr double temperature{20.0};              //!< T_eq,e = T_eq,i, C
        inline constexpr double exteriorHumidity{0.45};         //!< phi_a,e for t > 0
        inline constexpr double interiorHumidity{0.65};         //!< phi_a,i for t > 0
        inline constexpr double filmCoefficient{25.0};          //!< alpha_e, W/(m2 K)
        inline constexpr double vapourCoefficient{1.0e-3};      //!< beta_p, s/m, both sides
        inline constexpr double initialWater{84.7687};          //!< kg/m3 (RH = 95 %)
        inline constexpr double initialHumidity{0.95};

        //! The engine's resistance factor for the printed constant delta_p: its vapour flux
        //! is -(D_air(T) / mu) grad(c) with c = p M_w / (R T), so
        //! delta_p = D_air M_w / (mu R T) and mu = D_air,engine(T_ref) M_w / (R T_ref delta_p).
        [[nodiscard]] inline double engineResistanceFactor()
        {
            const double engineDiffusivityInAir{
              HygroThermFEM::vaporDiffusionCoefficientAtTemperature(temperature)};
            return engineDiffusivityInAir * molarMassWater
                   / (gasConstant * referenceTemperature * vapourPermeability);
        }

        //! Humidity samples for the isotherm table, dense across the benchmark's whole
        //! range (0,45 to 0,95) and closed at both ends. The analytic solution is linear
        //! in w, so the engine's humidity-based solve must see an isotherm whose secant
        //! capacity and tangent slope agree closely; 400 points hold both to a few 1e-4.
        [[nodiscard]] inline std::vector<FenestrationCommon::point>
          sorptionTable(const std::size_t points = 400u)
        {
            std::vector<FenestrationCommon::point> table{{0.0, 0.0}};
            table.reserve(points + 2u);
            for(std::size_t index = 1u; index < points; ++index)
            {
                const double humidity{static_cast<double>(index) / static_cast<double>(points)};
                table.emplace_back(humidity, waterContent(humidity));
            }
            table.emplace_back(1.0, 116.0);
            return table;
        }

        //! Material A as the engine consumes it. The constant D_w is the report's TOTAL
        //! moisture diffusivity, D_w = -K dP_suc/dw + delta_p p_s / xi ([B] p. 20); the
        //! vapour share delta_p p_s / xi is below 1e-14 m2/s over the benchmark's range,
        //! five orders under 6e-10, so the whole of it is handed to the liquid curve and
        //! the printed delta_p goes in as it stands.
        [[nodiscard]] inline HygroThermFEM::SolidMaterialParams material()
        {
            return {
              .name = "HAMSTAD BM2 material A",
              .thermalConductivityDry = thermalConductivity,
              .density = density,
              .porosity = 116.0 / densityWater,
              .heatCapacity = heatCapacity,
              .diffusionResistanceFactor = engineResistanceFactor(),
              .diffusionResistanceFactorMoistureDependent = {},
              .thermalConductivityMoistureDependent = {{0.0, thermalConductivity},
                                                       {1.0, thermalConductivity}},
              .moistureDependentMeasurementTemperature = 0,
              .thermalConductivityTemperatureDependent = {{0.0, thermalConductivity},
                                                          {100.0, thermalConductivity}},
              .temperatureDependentMeasurementHumidity = 0,
              .liquidTransportCurve = {{0.0, moistureDiffusivity}, {116.0, moistureDiffusivity}},
              .sorptionCurve = sorptionTable()};
        }
    }   // namespace Benchmark2

    /////////////////////////////////////////////////////////////////////////////
    /// Benchmark 5, "Capillary active inside insulation", [B] p. 45-49
    /////////////////////////////////////////////////////////////////////////////
    namespace Benchmark5
    {
        //! The retention scales are printed in 1/m against the suction height
        //! h = P_suc / (rho_w g), [B] p. 46; this turns one into the 1/Pa the modes carry.
        [[nodiscard]] inline constexpr double scalePerPascal(const double scalePerMetre)
        {
            return scalePerMetre / (densityWater * gravity);
        }

        //! One column of the parameter table on [B] p. 47.
        [[nodiscard]] inline Material brick()
        {
            return {.name = "HAMSTAD BM5 brick",
                    .freeSaturation = 373.5,
                    .retention = {{.weight = 0.46,
                                   .scale = scalePerPascal(0.47),
                                   .exponentN = 1.5,
                                   .exponentM = 1.0 - 1.0 / 1.5},
                                  {.weight = 0.54,
                                   .scale = scalePerPascal(0.2),
                                   .exponentN = 3.8,
                                   .exponentM = 1.0 - 1.0 / 3.8}},
                    .dryResistanceFactor = 7.5,
                    .deficitFactorP = 0.20,
                    .liquidConductivity = exponentialPolynomialInVolumetricWater(
                      {-36.484, 461.325, -5240.0, 2.907e4, -7.41e4, 6.997e4}),
                    .dryThermalConductivity = 0.682,
                    .conductivityMoistureSlope = 0.0 / densityWater,
                    .dryVolumetricHeatCapacity = 1600.0 * 1000.0};
        }

        [[nodiscard]] inline Material mortar()
        {
            return {.name = "HAMSTAD BM5 mortar",
                    .freeSaturation = 700.0,
                    .retention = {{.weight = 0.2,
                                   .scale = scalePerPascal(0.5),
                                   .exponentN = 1.5,
                                   .exponentM = 1.0 - 1.0 / 1.5},
                                  {.weight = 0.8,
                                   .scale = scalePerPascal(0.004),
                                   .exponentN = 3.8,
                                   .exponentM = 1.0 - 1.0 / 3.8}},
                    .dryResistanceFactor = 50.0,
                    .deficitFactorP = 0.20,
                    .liquidConductivity = exponentialPolynomialInVolumetricWater(
                      {-40.425, 83.319, -175.961, 123.863, 0.0, 0.0}),
                    .dryThermalConductivity = 0.6,
                    .conductivityMoistureSlope = 0.56 / densityWater,
                    .dryVolumetricHeatCapacity = 230.0 * 920.0};
        }

        [[nodiscard]] inline Material insulation()
        {
            return {.name = "HAMSTAD BM5 inside insulation",
                    .freeSaturation = 871.0,
                    .retention = {{.weight = 0.41,
                                   .scale = scalePerPascal(0.006),
                                   .exponentN = 2.5,
                                   .exponentM = 1.0 - 1.0 / 2.5},
                                  {.weight = 0.59,
                                   .scale = scalePerPascal(0.012),
                                   .exponentN = 2.4,
                                   .exponentM = 1.0 - 1.0 / 2.4}},
                    .dryResistanceFactor = 5.6,
                    .deficitFactorP = 0.20,
                    .liquidConductivity = exponentialPolynomialInVolumetricWater(
                      {-46.245, 294.506, -1439.0, 3249.0, -3370.0, 1305.0}),
                    .dryThermalConductivity = 0.06,
                    .conductivityMoistureSlope = 0.56 / densityWater,
                    .dryVolumetricHeatCapacity = 212.0 * 1000.0};
        }

        inline constexpr double brickThickness{0.365};        //!< m, [B] p. 47
        inline constexpr double mortarThickness{0.015};       //!< m
        inline constexpr double insulationThickness{0.040};   //!< m
        inline constexpr double interiorTemperature{20.0};    //!< C, [B] p. 48
        inline constexpr double interiorHumidity{0.60};
        inline constexpr double exteriorTemperature{0.0};   //!< C
        inline constexpr double exteriorHumidity{0.80};
        inline constexpr double exteriorFilmCoefficient{25.0};            //!< alpha_e,e W/(m2 K)
        inline constexpr double interiorFilmCoefficient{8.0};             //!< alpha_e,i
        inline constexpr double exteriorVapourCoefficient{1.838200e-7};   //!< beta_p,e s/m
        inline constexpr double interiorVapourCoefficient{5.8823e-8};     //!< beta_p,i s/m
        inline constexpr double initialTemperature{25.0};                 //!< C
        inline constexpr double initialHumidity{0.60};
        inline constexpr double simulationDays{60.0};
    }   // namespace Benchmark5
}   // namespace TestHelper::HAMSTAD
