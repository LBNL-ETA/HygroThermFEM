#pragma once

#include "IBCLine2D.hxx"
#include "Material.hxx"

namespace MoisThermFEM
{
    ////////////////////////////////////////////////////////
    /// ConvectionBC
    ////////////////////////////////////////////////////////

    //! \brief Convection boundary condition
    //!
    //! Simple convection boundary condition that is used in thermal module. It constant
    //! heat transfer coefficient and outside air temperature.
    class ConvectionBC : public IBCLinear2D
    {
    public:
        //! Constructor for convection boundary condition requires two nodes, heat transfer
        //! coefficient and air temperature.
        ConvectionBC(size_t index1,                    //!< Node 1 index
                     size_t index2,                    //!< Node 2 index
                     double t_ConvectionCoefficient,   //!< Film (convection) coefficient
                     double t_AirTemperature           //!< Outside air temperature
        );

        //! Function that calculates right hand side vector.
        std::vector<double> R_Vector() const override;

        //! Function that calculates matrix.
        FenestrationCommon::SquareMatrix H_Matrix() const override;

    protected:
        const double m_ConvectionCoefficient;
        const double m_AirTemperature;
    };

    ////////////////////////////////////////////////////////
    /// TemperatureBC
    ////////////////////////////////////////////////////////

    //! Constant temperature boundary condition
    //!
    //! TemperatureBC will be just special case of convection BC with huge value
    //! for film coefficients. It is used only in thermal module.
    class TemperatureBC : public ConvectionBC
    {
    public:
        //! Construction of temperature boundary condition with two nodes and identical temperatures
        //! in both of them
        TemperatureBC(size_t index1,              //!< Node 1 index
                      size_t index2,              //!< Node 2 index
                      double t_NodeTemperatures   //!< Temperature at both nodes
        );

        //! Construction of temperature boundary condition with two nodes and temperature at each of
        //! them
        TemperatureBC(size_t index1,    //!< Node 1 index
                      size_t index2,    //!< Node 2 index
                      double t_Temp1,   //!< Temperature at node 1
                      double t_Temp2    //!< Temperature at node 2
        );
    };

    ////////////////////////////////////////////////////////
    /// Flux BC
    ////////////////////////////////////////////////////////

    //! \brief Constant flux boundary condition.
    //!
    //! Boundary condition is used only in thermal module.
    class FluxBC : public IBCLinear2D
    {
    public:
        //! Construction of constant flux boundary condition with two nodes and constant flux.
        FluxBC(size_t index1,   //!< Node 1 index
               size_t index2,   //!< Node 2 index
               double t_Flux    //!< Constant flux at both nodes
        );

        //! Function that calculates right hand side vector.
        std::vector<double> R_Vector() const override;

        //! Function that calculates matrix.
        FenestrationCommon::SquareMatrix H_Matrix() const override;

    private:
        double m_Flux;
    };

    ///////////////////////////////////////////////////////
    /// BlackBodyRadiationBC
    ///////////////////////////////////////////////////////

    //! \brief Simple black body radiation boundary condition.
    //!
    //! Boundary condition is used only in thermal module.
    class BlackBodyRadiationBC : public IBCLinear2D
    {
    public:
        //! Black body radiation boundary condition with two nodes, emissivity and radiation
        //! temperature
        BlackBodyRadiationBC(
          size_t index1,                  //!< Node 1 index
          size_t index2,                  //!< Node 2 index
          double t_Emissivity,            //!< Boundary condition surface emissivity at both nodes.
          double t_RadiationTemperature   //!< Outside radiation temperature.
        );

        //! Function that calculates right hand side vector.
        std::vector<double> R_Vector() const override;

        //! Function that calculates matrix.
        FenestrationCommon::SquareMatrix H_Matrix() const override;

    private:
        //! Radiative convective coefficient that needs to be calculated based on current
        //! temperatures
        std::vector<double> HRadiative() const;

        double m_RadiationTemperature;
        double m_Emissivity;
    };


    /////////////////////////////////////////////////////
    /// MoistureBC
    /////////////////////////////////////////////////////

    //! \brief Moisture boundary condition with constant heat transfer coefficient, air humidity and
    //! air temperature.
    //!
    //! Boundary condition is used only in thermal and mass transfer module.
    class MoistureBC : public IBCLinear2D
    {
    public:
        //! Construction for moisture boundary condition.
        MoistureBC(size_t index1,   //!< Node 1 index
                   size_t index2,   //!< Node 2 index
                   const std::string &
                     materialName,   //!< Material name adjacent to boundary condition line.
                   double t_ConvectiveCoefficient,   //!< Heat transfer coefficient.
                   double t_AirHumidity,             //!< Outside air humidity.
                   double t_AirTemperature           //!< Outside air temperature.
        );

        //! Function that calculates right hand side vector.
        std::vector<double> R_Vector() const override;

        //! Function that calculates matrix.
        FenestrationCommon::SquareMatrix H_Matrix() const override;

    protected:
        double m_ConvectiveCoefficient;
        double m_AirHumidity;
        double m_AirTemperature;
        const Material & m_Material;
    };

}   // namespace MoisThermFEM
