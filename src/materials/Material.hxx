#pragma once

#include <vector>
#include <memory>
#include <string>
#include <optional>

#include "Functions.hxx"
#include "KeffCavity.hxx"

namespace HygroThermFEM
{
    //! \brief Enumerates different state of water content.
    //!
    //! Water content can be obtained in three different states or simply in total value of water
    //! at given point.
    enum class WaterContent
    {
        Water,    //!< Total water content.
        Liquid,   //!< Water content in liquid state.
        Vapor,    //!< Water content in vapor state.
        Ice       //!< Water content in frozen state.
    };

    class Water
    {
    public:
        Water(double water = 0, double liquid = 0, double vapor = 0, double ice = 0);

        Water & operator*(const double & other);
        Water & operator+=(const Water & other);
        Water & operator/(const double & other);

        [[nodiscard]] double content(WaterContent content) const;

    private:
        std::map<WaterContent, double> m_Content;
    };

    //! \brief Standard used in thermal calculations of air pockets.
    enum class CavityStandard
    {
        ISO15099,
        CEN
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // IMaterial
    ///////////////////////////////////////////////////////////////////////////////////////////////

    //! \brief Interface class for all materials that will be used in analysis.
    class IMaterial
    {
    public:
        virtual ~IMaterial() = default;

        IMaterial(
          std::string name,
          double thermalConductivityDry,
          double density,
          double porosity,
          double heatCapacity,
          double diffusionResistanceFactor,
          const std::vector<FenestrationCommon::point> & thermalConductivityMoistureDependent,
          double moistureDependentMeasurementTemperature,
          const std::vector<FenestrationCommon::point> & thermalConductivityTemperatureDependent,
          double temperatureDependentMeasurementHumidity,
          const std::vector<FenestrationCommon::point> & liquidTransportationCurve,
          const std::vector<FenestrationCommon::point> & sorptionCurve,
          double emissivity,
          bool isLinear = true);

        IMaterial(std::string name, bool isLinear = true);

        //! Material's name.
        [[nodiscard]] std::string name() const;

        //! Thermal conductivity of dry material
        [[nodiscard]] virtual double thermalConductivityDry() const final;

        //! Sets material thermal conductivity to given value
        //!
        //! \param thermalConductivity Value to which thermal conductivity will be set to
        virtual void setThermalConductivity(double thermalConductivity) final;

        //! \brief Checks if value has been set
        [[nodiscard]] virtual bool hasThermalConductivityDry() const final;

        //! Material's density.
        [[nodiscard]] virtual double density() const final;

        //! Sets material density to given value
        //!
        //! \param density Value to which density will be set to
        virtual void setDensity(double density) final;

        //! \brief Checks if value has been set
        [[nodiscard]] virtual bool hasDensity() const final;

        //! Material's specific heat capacity.
        [[nodiscard]] virtual double heatCapacity() const final;

        //! Sets material heat capacity to given value
        //!
        //! \param heatCapacity Value to which heat capacity will be set to
        virtual void setHeatCapacity(double heatCapacity) final;

        //! \brief Checks if value has been set
        [[nodiscard]] virtual bool hasHeatCapacity() const final;

        //! Material's porosity.
        [[nodiscard]] virtual double porosity() const final;

        //! Sets material porosity to given value
        //!
        //! \param porosity Value to which porosity will be set to
        virtual void setPorosity(double porosity) final;

        //! \brief Checks if value has been set
        [[nodiscard]] virtual bool hasPorosity() const final;

        //! Returns material emissivity.
        [[nodiscard]] virtual double emissivity() const final;

        //! Sets material emissivity to given value
        //!
        //! \param emissivity Value to which emissivity will be set to
        virtual void setEmissivity(double emissivity) final;

        //! \brief Checks if value has been set
        [[nodiscard]] virtual bool hasEmissivity() const final;

        //! Material's diffusion resistance factor.
        [[nodiscard]] virtual double diffusionResistanceFactor() const final;

        //! Sets material diffusion resistance factor to given value
        //!
        //! \param diffusionResistanceFactor Value to which diffusion resistance factor will be set
        //! to
        virtual void setDiffusionResistanceFactor(double diffusionResistanceFactor) final;

        //! \brief Checks if value has been set
        [[nodiscard]] virtual bool hasDiffusionResistanceFactor() const final;

        //! Material can require nonlinear iterations.
        [[nodiscard]] virtual bool isLinear() const final;

        //! Thermal conductivity table (x-water content [kg/m3], y-thermal conductivity[W/(mK)])
        [[nodiscard]] virtual TabularFunction2D thermalConductivityMoistureAndTemperatureDependent() const final;

        //! Sets thermal condctivity dependence on temperature and moisture
        //!
        //! \param thermalConductivityMoistureDependent table of thermal conductivity values in
        //! relation to moisture content
        //! \param moistureDependentMeasuredTemperature temperature at
        //! which moisture dependent thermal conductivity is measured
        //! \param thermalConductivityTemperatureDependent table of thermal conductivity values
        //! in relation to temperature
        //! \param temperatureDependentMeasurementHumidity moisture content
        //! at which thermal conductivity table is measured
        virtual void setThermalConductivityMoistureAndTemperatureDependent(
          const std::vector<FenestrationCommon::point> & thermalConductivityMoistureDependent,
          double moistureDependentMeasuredTemperature,
          const std::vector<FenestrationCommon::point> & thermalConductivityTemperatureDependent,
          double temperatureDependentMeasurementHumidity) final;

        //! \brief Checks if table has been set
        [[nodiscard]] virtual bool hasThermalConductivityMoistureAndTemperatureDependent() const final;

        //! \brief Liquid transportation curve of the material.
        //!
        //! Liquid transportation coefficient shows how
        //! much of water can be distributed through the material with certain water content
        //! (x-water content [kg/m3], y-water flow [m2/s]
        [[nodiscard]] virtual const std::vector<FenestrationCommon::point> &
          liquidTransportationCurve() const final;

        //! \brief Sets material liquid transportation curve
        //!
        //! \param liquidTransportationCurve values at which liquid transportation curve will be set
        //! to
        virtual void setLiquidTransportationCurve(
          const std::vector<FenestrationCommon::point> & liquidTransportationCurve) final;

        //! \brief Checks if curve has been set
        [[nodiscard]] virtual bool hasLiquidTransportationCurve() const final;

        //! \brief Material's sorption curve.
        //!
        //! Sorption curve or moisture storage function show how
        //! much of water content is contained in the material at certain relative humidity
        //! (x-relative humidity [-], y-water content [kg/m3])
        [[nodiscard]] virtual const std::vector<FenestrationCommon::point> & sorptionCurve() const final;

        //! \brief Sets material's sorption curve
        //!
        //! \param sorptionCurve values for sorption curve
        virtual void
          setSorptionCurve(const std::vector<FenestrationCommon::point> & sorptionCurve) final;

        //! \brief Checks if curve has been set
        [[nodiscard]] virtual bool hasSorptionCurve() const final;

        //! \brief Water content in given node.
        //!
        //! \param node Node for which water content is required.
        //! \param waterContent Water content property (total, liquid, vapor or ice).
        //! \return Value of water content.
        [[nodiscard]] virtual Water waterContent(const INode2D & node) const = 0;

        // Materials are stored in set that is used by other classes and that is why we need these
        // operators
        friend bool operator<(const IMaterial & lhs, const IMaterial & rhs);
        friend bool operator>(const IMaterial & lhs, const IMaterial & rhs);
        friend bool operator<=(const IMaterial & lhs, const IMaterial & rhs);
        friend bool operator>=(const IMaterial & lhs, const IMaterial & rhs);

    protected:

        //! Saturated vapor content calculations at given node. It is necessary
        //! for water content calculations.
        static double saturationConcentration(const INode2D & node);

        std::string m_Name;
        std::optional<double> m_ThermalConductivityDry{};
        std::optional<double> m_Density{};
        std::optional<double> m_Porosity{};
        std::optional<double> m_SpecificHeatCapacity{};
        std::optional<double> m_DiffusionResistanceFactor{};

        //! Thermal conductivity table is (x-water content [kg/m3], temperature[Celsius], y-thermal
        //! conductivity[W/(mK)])
        std::unique_ptr<TabularFunction2D> m_ThermalConductivity2DTable;

        //! Liquid transportation coefficient is function of water content. It shows how much of
        //! water will be transferred through material in relation to water content (x-water content
        //! [kg/m3], y-liquid transportation coefficient [m2/s]
        std::unique_ptr<TabularFunction1D> m_LiquidTransportCoefficient;

        //! Sorption curve shows how much of water content will be in relation to relative humidity
        //! (x-relative humidity [between zero to one], y-water content [kg/m3]
        std::unique_ptr<TabularFunction1D> m_SorptionCurve;

        std::optional<double> m_Emissivity{};
        const bool m_Linear;
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // IGas
    ///////////////////////////////////////////////////////////////////////////////////////////////

    //! Interface for gases
    class IGas : public IMaterial
    {
    public:
        IGas(std::string name,
             CavityStandard cavityStandard, Gases::CGas gas = Gases::CGas());

        //! \brief Some materials will require update of thermal conductivity within iterations.
        //! This virtual function requires update in every material type.
        //!
        //! \param thermalConductivity New value for thermal conductivity.
        virtual void updateThermalConductivity(double thermalConductivity) = 0;

        //! \brief Gases will require update on diffusion resistance factor on every iteration.
        //!
        //! \param diffusionResistanceFactor New value for diffusion resistance factor.
        virtual void updateDiffusionResistanceFactor(double diffusionResistanceFactor) = 0;

        //! \brief Update sorption curve with new data
        //!
        //! \param saturationTemperature Temperature at which saturation will occur
        virtual void updateSorptionCurve(double saturationTemperature) = 0;

        [[nodiscard]] CavityStandard standard() const;

        [[nodiscard]] Gases::CGas & getGas();

    private:
        CavityStandard m_CavityStandard;

        Gases::CGas m_Gas;

        static const double DefaultThermalConductivityDry;
        static const double DefaultDensity;
        static const double DefaultPorosity;
        static const double DefaultHeatCapacity;
        const static double DefaultDiffusionResistanceFactor;
        const static std::vector<FenestrationCommon::point> DefaultSorptionCurve;
        const static std::vector<FenestrationCommon::point> DefaultThermalConductivityMoistureDependent;
        static const double DefaultMoistureDependentMeasurementTemperature;
        const static std::vector<FenestrationCommon::point> DefaultThermalConductivityTemperatureDependent;
        static const double DefaultTemperatureDependentMeasurementMoisture;
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // SolidMaterialParams
    ///////////////////////////////////////////////////////////////////////////////////////////////

    //! \brief Parameters for creating a SolidMaterial with designated initializers (C++20)
    struct SolidMaterialParams
    {
        std::string name;
        double thermalConductivityDry = 0.0;
        double density = 0.0;
        double porosity = 0.0;
        double heatCapacity = 0.0;
        double diffusionResistanceFactor = 0.0;
        std::vector<FenestrationCommon::point> thermalConductivityMoistureDependent = {};
        double moistureDependentMeasurementTemperature = 0.0;
        std::vector<FenestrationCommon::point> thermalConductivityTemperatureDependent = {};
        double temperatureDependentMeasurementHumidity = 0.0;
        std::vector<FenestrationCommon::point> liquidTransportCurve = {};
        std::vector<FenestrationCommon::point> sorptionCurve = {};
        double emissivity = 0.9;
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // SolidMaterial
    ///////////////////////////////////////////////////////////////////////////////////////////////

    //! \brief Keeps necessary material data and performs some minor water content calculations.
    //!
    //! SolidMaterial data necessary for mass and thermal transfer are held in this class. Mainly
    //! used as storage.
    class SolidMaterial : public IMaterial
    {
    public:
        SolidMaterial() = delete;

        //! \brief SolidMaterial construction is done through singleton class
        //!
        //! \param name SolidMaterial name
        //! \param density Material density
        //! \param porosity Material porosity
        //! \param heatCapacity Specific heat capacity of dry material
        //! \param diffusionResistanceFactor Diffuse resistance factor
        //! \param thermalConductivityMoistureDependent Moisture dependent thermal conductivity
        //! \param moistureDependentMeasurementTemperature Temperature at which moisture dependent
        //! thermal conductivity is measured
        //! \param thermalConductivityTemperatureDependent Temperature dependent thermal
        //! conductivity
        //! \param temperatureDependentMeasurementHumidity Humidity at which
        //! temperature dependent thermal conductivity is measured
        //! \param liquidTransportCurve Liquid transportation curve. Relationship between relative
        //! humidity and ability of material to transport water.
        //! \param sorptionCurve Moisture storage function.
        //! Relationship between relative humidity and water content.
        //! \param emissivity SolidMaterial emissivity
        SolidMaterial(
          std::string name,
          double thermalConductivityDry,
          double density,
          double porosity,
          double heatCapacity,
          double diffusionResistanceFactor,
          const std::vector<FenestrationCommon::point> & thermalConductivityMoistureDependent,
          double moistureDependentMeasurementTemperature,
          const std::vector<FenestrationCommon::point> & thermalConductivityTemperatureDependent,
          double temperatureDependentMeasurementHumidity,
          const std::vector<FenestrationCommon::point> & liquidTransportCurve,
          const std::vector<FenestrationCommon::point> & sorptionCurve,
          double emissivity = 0.9);

        explicit SolidMaterial(std::string name);

        //! \brief Construct SolidMaterial from params struct (C++20 designated initializers)
        explicit SolidMaterial(SolidMaterialParams params);

        //! \brief Water content in given node.
        //!
        //! \param node Node for which water content is required.
        //! \param waterContent Water content property (total, liquid, vapor or ice).
        //! \return Value of water content.
        [[nodiscard]] Water waterContent(const INode2D & node) const override;

    private:
        //! Returns total water content in given node.
        [[nodiscard]] double totalWaterContent(const INode2D & node) const;

        //! Returns vapor content in given node.
        [[nodiscard]] double vaporContent(const INode2D & node) const;

        //! Returns liquid content in given node.
        [[nodiscard]] double liquidWaterContent(const INode2D & node) const;

        //! Returns ice content in given node.
        [[nodiscard]] double iceContent(const INode2D & node) const;

        //! Calculates amount of pores filled with liquid. Necessary for water content calculations.
        [[nodiscard]] double liquidPorosity(const INode2D & node) const;

        //! Calculates amount of pores filled with air. Necessary for water content calculations.
        [[nodiscard]] double airPorosity(const INode2D & node) const;
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Gas
    ///////////////////////////////////////////////////////////////////////////////////////////////

    class Gas : public IGas
    {
        friend class Materials;

    public:
        Gas() = delete;

        explicit Gas(const std::string & name, CavityStandard cavityStandard = CavityStandard::ISO15099, Gases::CGas gas = Gases::CGas());

        //! Water content for given node
        //!< param node Node for which water content is required.
        [[nodiscard]] Water waterContent(const INode2D & node) const override;

        void updateThermalConductivity(double thermalConductivity) override;

        void updateDiffusionResistanceFactor(double diffusionResistanceFactor) override;

        void updateSorptionCurve(double saturationTemperature) override;        
    };
}   // namespace HygroThermFEM
