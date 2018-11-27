#pragma once

#include <vector>
#include <memory>
#include <string>

#include "Functions.hxx"

namespace MoisThermFEM
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

        //! Material's thermal conductivity.
        double thermalConductivity() const;

        //! Material's diffusion resistance factor.
        double diffusionResistanceFactor() const;
        std::vector<std::pair<double, double>> liquidTransportationCurve() const;

        //! Water content for given node
        double waterContent(const Node2D & node,   //!< Node for which water content is required.
                            WaterContent waterContent   //!< Water content property.
                            ) const;

        //! Material's sorption curve
        std::vector<std::pair<double, double>> sorptionCurve() const;

        /// Materials will be stored in set which require operator >.

        friend bool operator<(const Material & lhs, const Material & rhs);
        friend bool operator>(const Material & lhs, const Material & rhs);
        friend bool operator<=(const Material & lhs, const Material & rhs);
        friend bool operator>=(const Material & lhs, const Material & rhs);

    private:
        //! Returns total water content in given node.
        double waterContent(const Node2D & node) const;

        //! Returns vapor content in given node.
        double vaporContent(const Node2D & node) const;

        //! Returns liquid content in given node.
        double liquidWaterContent(const Node2D & node) const;

        //! Returns ice content in given node.
        double iceContent(const Node2D & node) const;

        /// Create material by using MaterialPool.

        //! \brief Material construction is done through singleton class
        Material(const std::string & Name,
                 double Density,
                 double Porosity,
                 double HeatCapacity,
                 double ThermalConductivity,
                 double DiffusionResistanceFactor,
                 const std::vector<std::pair<double, double>> & LiquidTransportCurve,
                 const std::vector<std::pair<double, double>> & SorptionCurve);

        std::string m_Name;
        double m_Density;
        double m_Porosity;
        double m_HeatCapacity;
        double m_ThermalConductivity;
        double m_DiffusionResistanceFactor;
        std::unique_ptr<MoisThermFEM::TabularFunction> m_LiquidTransportCoefficient;
        std::unique_ptr<MoisThermFEM::TabularFunction> m_SorptionCurve;

        double saturatedVaporContent(const Node2D & node) const;
        double liquidPorosity(const Node2D & node) const;
        double airPorosity(const Node2D & node) const;
    };

}   // namespace MoisThermFEM
