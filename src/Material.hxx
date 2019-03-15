#pragma once

#include <vector>
#include <memory>
#include <string>

#include "Functions.hxx"

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

    //! \brief Used to describe different type of materials.
    //!
    //! Depending on material type, engine will perform different algorithms on how to calculate
    //! equivalent material properties.
    enum class MaterialType
    {
        Solid,
        Gas
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // IMaterial
    ///////////////////////////////////////////////////////////////////////////////////////////////

    //! \brief Interface class for all materials that will be used in analysis.
    class IMaterial
    {
    public:
        virtual ~IMaterial() = default;

        IMaterial(std::string cs,
                  double density,
                  double porosity,
                  double heatCapacity,
                  double diffusionResistanceFactor,
                  const std::vector<std::pair<double, double>> & thermalConductivity,
                  const std::vector<std::pair<double, double>> & liquidTransportationCurve,
                  const std::vector<std::pair<double, double>> & sorptionCurve,
                  double emissivity,
                  MaterialType material);

        //! Material's name.
        std::string name() const;

        //! Material's density.
        double density() const;

        //! Material's specific heat capacity.
        double heatCapacity() const;

        //! Material's porosity.
        double porosity() const;

        //! Returns type of material.
        MaterialType materialType() const;

        //! Returns material emissivity.
        double emissivity() const;

        //! Material's diffusion resistance factor.
        double diffusionResistanceFactor() const;

        //! Thermal conductivity table (x-water content [kg/m3], y-thermal conductivity[W/(mK)])
        const std::vector<std::pair<double, double>> & thermalConductivity() const;

        //! \brief Liquid transportation curve of the material.
        //!
        //! Liquid transportation coefficient shows how
        //! much of water can be distributed through the material with certain water content
        //! (x-water content [kg/m3], y-water flow [m2/s]
        const std::vector<std::pair<double, double>> & liquidTransportationCurve() const;

        //! \brief Material's sorption curve.
        //!
        //! Sorption curve or moisture storage function show how
        //! much of water content is contained in the material at certain relative humidity
        //! (x-relative humidity [-], y-water content [kg/m3])
        const std::vector<std::pair<double, double>> & sorptionCurve() const;

        //! \brief Water content in given node.
        //!
        //! \param node Node for which water content is required.
        //! \param waterContent Water content property (total, liquid, vapor or ice).
        //! \return Value of water content.
        virtual double waterContent(const INode2D & node, WaterContent waterContent) const = 0;

        //! \brief Some materials will require update of thermal conductivity within iterations.
        //! This virtual function requires update in every material type.
        //!
        //! \param thermalConductivity New value for thermal conductivity.
        virtual void updateThermalConductivity(double thermalConductivity) = 0;

        friend bool operator<(const IMaterial & lhs, const IMaterial & rhs);
        friend bool operator>(const IMaterial & lhs, const IMaterial & rhs);
        friend bool operator<=(const IMaterial & lhs, const IMaterial & rhs);
        friend bool operator>=(const IMaterial & lhs, const IMaterial & rhs);

    protected:
        std::string m_Name;
        double m_Density;
        double m_Porosity;
        double m_HeatCapacity;
        double m_DiffusionResistanceFactor;

        //! Thermal conductivity table is (x-water content [kg/m3], y-thermal conductivity[W/(mK)])
        std::unique_ptr<TabularFunction> m_ThermalConductivity;

        //! Liquid transportation coefficient is function of water content. It shows how much of
        //! water will be transferred through material in relation to water content (x-water content
        //! [kg/m3], y-liquid transportation coefficient [m2/s]
        std::unique_ptr<TabularFunction> m_LiquidTransportCoefficient;

        //! Sorption curve shows how much of water content will be in relation to relative humidity
        //! (x-relative humidity [between zero to one], y-water content [kg/m3]
        std::unique_ptr<TabularFunction> m_SorptionCurve;

        double m_Emissivity;
        MaterialType m_MaterialType;
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
        friend class MaterialPool;

    public:
        SolidMaterial() = delete;

        //! \brief Water content in given node.
        //!
        //! \param node Node for which water content is required.
        //! \param waterContent Water content property (total, liquid, vapor or ice).
        //! \return Value of water content.
        double waterContent(const INode2D & node, WaterContent waterContent) const override;

        void updateThermalConductivity(double thermalConductivity) override;

    private:
        //! Returns total water content in given node.
        double waterContent(const INode2D & node) const;

        //! Returns vapor content in given node.
        double vaporContent(const INode2D & node) const;

        //! Returns liquid content in given node.
        double liquidWaterContent(const INode2D & node) const;

        //! Returns ice content in given node.
        double iceContent(const INode2D & node) const;

        //! \brief SolidMaterial construction is done through singleton class
        SolidMaterial(
          const std::string & name,           //!< SolidMaterial name
          double density,                     //!< Density of dry material
          double porosity,                    //!< SolidMaterial porosity
          double heatCapacity,                //!< Specific heat capacity of dry material
          double diffusionResistanceFactor,   //!< Diffuse resistance factor
          const std::vector<std::pair<double, double>> &
            thermalConductivity,   //!< SolidMaterial conductivity of dry material where
                                   //!< conductivity depends on water content
          const std::vector<std::pair<double, double>> &
            liquidTransportCurve,   //!< Liquid transportation curve. Relationship between relative
                                    //!< humidity and ability of material to transport water.
          const std::vector<std::pair<double, double>> &
            sorptionCurve,   //!< Moisture storage function. Relationship between relative humidity
                             //!< and water content.
          double emissivity = 0.9   //!< SolidMaterial emissivity
        );

        //! Saturated vapor content calculations at given node. It is necessary
        //! for water content calculations.
        static double saturationConcentration(const INode2D & node);

        //! Calculates amount of pores filled with liquid. Necessary for water content calculations.
        double liquidPorosity(const INode2D & node) const;

        //! Calculates amount of pores filled with air. Necessary for water content calculations.
        double airPorosity(const INode2D & node) const;
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Gas
    ///////////////////////////////////////////////////////////////////////////////////////////////

    class Gas : public IMaterial
    {
        friend class MaterialPool;

    public:
        Gas() = delete;

        Gas(const std::string & name);

        //! Water content for given node
        double waterContent(const INode2D & node,   //!< Node for which water content is required.
                            WaterContent waterContent   //!< Water content property (total water,
                                                        //!< liquid, vapor or ice).
                            ) const override;

        void updateThermalConductivity(double thermalConductivity) override;
    };
}   // namespace HygroThermFEM
