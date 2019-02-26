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
        FrameCavity_ISO15099
    };

    //! \brief Keeps necessary material data
    //!
    //! Material data necessary for mass and thermal transfer are held in this class. Mainly used as
    //! storage. No heavy calculations are performed by this class.
    class Material
    {
        friend class MaterialPool;

    public:
        Material() = delete;

        //! Material's name.
        std::string name() const;

        //! Material's density.
        double density() const;

        //! Material's specific heat capacity.
        double heatCapacity() const;

        //! Material's porosity.
        double porosity() const;

        //! Material's diffusion resistance factor.
        double diffusionResistanceFactor() const;

        //! Thermal conductivity table is (x-water content [kg/m3], y-thermal conductivity[W/(mK)])
        const std::vector<std::pair<double, double>> & thermalConductivity() const;

        //! Liquid transportation curve of the material. Liquid transportation coefficient shows how
        //! much of water can be distributed through the material with certain water content
        //! (x-water content [kg/m3], y-water flow [m2/s]
        const std::vector<std::pair<double, double>> & liquidTransportationCurve() const;

        //! Water content for given node
        double waterContent(const INode2D & node,   //!< Node for which water content is required.
                            WaterContent waterContent   //!< Water content property (total water,
                                                        //!< liquid, vapor or ice).
                            ) const;

        //! Material's sorption curve. Sorption curve or moisture storage function show how much of
        //! water content is contained in the material at certain relative humidity (x-relative
        //! humidity [-], y-water content [kg/m3])
        const std::vector<std::pair<double, double>> & sorptionCurve() const;

        /// Materials will be stored in set which require operator >.

        friend bool operator<(const Material & lhs, const Material & rhs);
        friend bool operator>(const Material & lhs, const Material & rhs);
        friend bool operator<=(const Material & lhs, const Material & rhs);
        friend bool operator>=(const Material & lhs, const Material & rhs);

    private:
        //! Returns total water content in given node.
        double waterContent(const INode2D & node) const;

        //! Returns vapor content in given node.
        double vaporContent(const INode2D & node) const;

        //! Returns liquid content in given node.
        double liquidWaterContent(const INode2D & node) const;

        //! Returns ice content in given node.
        double iceContent(const INode2D & node) const;

        /// Create material by using MaterialPool.

        //! \brief Material construction is done through singleton class
        Material(
          const std::string & Name,           //!< Material name
          double Density,                     //!< Density of dry material
          double Porosity,                    //!< Material porosity
          double HeatCapacity,                //!< Specific heat capacity of dry material
          double DiffusionResistanceFactor,   //!< Diffuse resistance factor
          const std::vector<std::pair<double, double>> &
            ThermalConductivity,   //!< Material conductivity of dry material where conductivity
                                   //!< depends on water content
          const std::vector<std::pair<double, double>> &
            LiquidTransportCurve,   //!< Liquid transportation curve. Relationship between relative
                                    //!< humidity and ability of material to transport water.
          const std::vector<std::pair<double, double>> &
            SorptionCurve,    //!< Moisture storage function. Relationship between relative humidity
                              //!< and water content.
          double emissivity = 0.9, //!< Material emissivity
          MaterialType materialType = MaterialType::Solid);

        std::string m_Name;
        double m_Density;
        double m_Porosity;
        double m_HeatCapacity;
        double m_DiffusionResistanceFactor;

        // Thermal conductivity table is (x-water content [kg/m3], y-thermal conductivity[W/(mK)])
        std::unique_ptr<HygroThermFEM::TabularFunction> m_ThermalConductivity;

        // Liquid transportation coefficient is function of water content. It shows how much of
        // water will be transferred through material in relation to water content (x-water content
        // [kg/m3], y-liquid transportation coefficient [m2/s]
        std::unique_ptr<HygroThermFEM::TabularFunction> m_LiquidTransportCoefficient;

        // Sorption curve shows how much of water content will be in relation to relative humidity
        // (x-relative humidity [between zero to one], y-water content [kg/m3]
        std::unique_ptr<HygroThermFEM::TabularFunction> m_SorptionCurve;

        double m_Emissivity;
        MaterialType m_MaterialType;

        //! Saturated vapor content calculations at given node. It is necessary
        //! for water content calculations.
        static double saturationConcentration(const INode2D &node);

        //! Calculates amount of pores filled with liquid. Necessary for water content calculations.
        double liquidPorosity(const INode2D & node) const;

        //! Calculates amount of pores filled with air. Necessary for water content calculations.
        double airPorosity(const INode2D & node) const;
    };

}   // namespace HygroThermFEM
