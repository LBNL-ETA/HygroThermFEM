#pragma once

#include "IBCLine2D.hxx"
#include "Material.hxx"
#include "BoundaryConditionCoefficients.hxx"

namespace HygroThermFEM
{
    class IConvectiveCoefficient;

    ////////////////////////////////////////////////////////
    /// IConvectionBC
    ////////////////////////////////////////////////////////

    //! \brief Interface boundary condition class that is used to create various convection boundary
    //! condition models
    //!
    //! Simple convection boundary condition that is used in thermal module.
    class IConvectionBC : public IBCLinear2D
    {
    public:
        //! Constructor for convection boundary condition that is common between models
        IConvectionBC(size_t index1,
                      size_t index2,
                      double t_AirTemperature,
                      std::unique_ptr<IConvectiveCoefficient>,
                      double t_AirHumidity = 0,
                      bool t_SimulateMoisture = true);

        //! Function that calculates right hand side vector.
        [[nodiscard]] std::vector<double> R_Vector() const override;

        //! Function that calculates H matrix.
        [[nodiscard]] SquareMatrix H_Matrix() const override;

    protected:
        const double m_AirTemperature;
        std::unique_ptr<IConvectiveCoefficient> m_ConvectiveCoeffCalc;
        const double m_AirHumidity;

        //! Need to include energy from vapor flux energy if moisture is present.
        bool m_SimulateVaporFluxEnergy;
    };

    ////////////////////////////////////////////////////////
    /// ConstantConvectionBC
    ////////////////////////////////////////////////////////
    class ConstantConvectionBC : public IConvectionBC
    {
    public:
        //! \brief Creation of constant convection heat transfer coefficient boundary conditions.
        //!
        //! @param index1 node1 index
        //! @param index2 node2 index
        //! @param fixedBCHCCoefficients Structure that holds necessary coefficients in other to
        //! create constant film coefficient boundary condition
        //! @param simulateVaporFluxEnergy Need to be true if we want to include vapor flux energy
        ConstantConvectionBC(size_t index1,
                             size_t index2,
                             const FixedBCHCCoefficients & fixedBCHCCoefficients,
                             bool simulateVaporFluxEnergy = true);
    };

    ////////////////////////////////////////////////////////
    /// TARPConvectionBC
    ////////////////////////////////////////////////////////
    class TARPConvectionBC : public IConvectionBC
    {
    public:
        //! \brief Convection boundary condition with TARP convection algorithm
        //!
        //! @param index1 node1 index
        //! @param index2 node2 index
        //! @param varHCCoeff Necessary coefficients that will be used to create TARP heat flow
        //! convection algorithm
        //! @param surfaceTilt Surface tilt at the boundary
        //! @param simulateVaporFluxEnergy Indicates whether or not moisture will be calculated at
        //! the boundary
        TARPConvectionBC(size_t index1,
                         size_t index2,
                         const VariableBCTARPHCCoefficients & varHCCoeff,
                         double surfaceTilt = 90,
                         bool simulateVaporFluxEnergy = true);
    };

    ////////////////////////////////////////////////////////
    /// TemperatureBC
    ////////////////////////////////////////////////////////

    //! Constant temperature boundary condition
    //!
    //! TemperatureBC will be just special case of convection BC with huge convectiveCoefficients
    //! for film coefficients. It is used only in thermal module.
    class TemperatureBC : public ConstantConvectionBC
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
        SquareMatrix H_Matrix() const override;

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
          double t_RadiationTemperature   //!< Outside radiation temperature in celsius.
        );

        //! Function that calculates right hand side vector.
        std::vector<double> R_Vector() const override;

        //! Function that calculates matrix.
        SquareMatrix H_Matrix() const override;

    private:
        //! Radiative convective coefficient that needs to be calculated based on current
        //! temperatures
        std::vector<double> HRadiative() const;

        double m_RadiationTemperature;
        double m_Emissivity;
    };

    ///////////////////////////////////////////////////////
    /// SimplifiedRadiationBC
    ///////////////////////////////////////////////////////
    class LinearizedRadiationBC : public IBCLinear2D
    {
    public:
        LinearizedRadiationBC(size_t index1,
                              size_t index2,
                              const LinearizedRadiationBCCoefficients & linearRadBC);

        std::vector<double> R_Vector() const override;

        SquareMatrix H_Matrix() const override;

    private:
        double m_RadiationCoefficient;
        double m_RadiationTemperature;
    };

    /////////////////////////////////////////////////////
    /// MoistureBC
    /////////////////////////////////////////////////////

    //! \brief Base class for moisture boundary conditions. The only difference in inherited classes
    //! is how convection heat transfer coefficient is calculated.
    //!
    //! This class handles all matrices and vectors generation that is common for moisture bcs
    class IMoistureBC : public IBCLinear2D
    {
    public:
        //! Construction for moisture boundary condition.
        IMoistureBC(size_t index1,
                    size_t index2,
                    const std::string & materialName,
                    double t_AirHumidity,
                    double t_AirTemperature,
                    std::unique_ptr<IConvectiveCoefficient> model);

        //! Function that calculates right hand side vector.
        [[nodiscard]] std::vector<double> R_Vector() const override;

        //! Function that calculates matrix.
        [[nodiscard]] SquareMatrix H_Matrix() const override;

    protected:
        double m_AirHumidity;
        double m_AirTemperature;
        const IMaterial & m_Material;

        //! \brief Convective coefficients calculations can be performed with different algorithms
        std::unique_ptr<IConvectiveCoefficient> m_ConvectiveCoeffCalc;
    };

    /////////////////////////////////////////////////////
    /// MoistureBCTARPHc
    /////////////////////////////////////////////////////

    //! \brief Moisture boundary condition that calculates convective coefficient.
    class MoistureBCTARPHc : public IMoistureBC
    {
    public:
        //! \brief Construction of boundary condition with TARP algorithm for convective heat
        //! transfer calculations
        //!
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param materialName Material name for which boundary is associated with
        //! @param varHCCoeff Coefficients that are necessary for TARP heat transfer coefficient
        //! calculations. Structure contain only coefficients that are variable through every
        //! timestep.
        //! @param surafceTilt Surface tilt at the boundary. [degrees]
        MoistureBCTARPHc(size_t index1,
                         size_t index2,
                         const std::string & materialName,
                         const VariableBCTARPHCCoefficients & varHCCoeff,
                         double surfaceTilt = 90);
    };

    /////////////////////////////////////////////////////
    /// MoistureBCFixedHc
    /////////////////////////////////////////////////////

    //! \brief Moisture boundary condition with fixed convective coefficient
    class MoistureBCFixedHc : public IMoistureBC
    {
    public:
        MoistureBCFixedHc(size_t index1,
                          size_t index2,
                          const std::string & materialName,
                          const FixedBCHCCoefficients & fixedBchcCoefficients);
    };

}   // namespace HygroThermFEM
