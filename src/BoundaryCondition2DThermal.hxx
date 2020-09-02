#pragma once

#include "BoundaryCondition2D.hxx"

namespace HygroThermFEM
{
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
    /// ThermalTARPConvectionBC
    ////////////////////////////////////////////////////////
    class ThermalTARPConvectionBC : public IConvectionBC
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
        ThermalTARPConvectionBC(size_t index1,
                                size_t index2,
                                const TARPCoefficients & varHCCoeff,
                                double surfaceTilt = 90,
                                bool simulateVaporFluxEnergy = true);
    };

    ////////////////////////////////////////////////////////
    /// ASHRAEInsideConvectionBC
    ////////////////////////////////////////////////////////

    //! \brief Convection boundary condition with ASHRAE inside convection algorithm
    class ASHRAEInsideConvectionBC : public IConvectionBC
    {
    public:
        //! \brief Construction of ASHRAE inside convection algorithm
        //!
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param coeff Timestep coefficients that are used in ASHRAE Inside convection algorithm
        //! @param surfaceHeight Surface height for which convection coefficient has been calculated
        //! [meters]
        //! @param surfaceTilt Surface tilt at the boundary [degrees]
        //! @param simulateVaporFluxEnergy Indicates whether energy from vapor flow should be
        //! included in calculations.
        ASHRAEInsideConvectionBC(size_t index1,
                                 size_t index2,
                                 const ASHRAEInsideCoefficients & coeff,
                                 double surfaceHeight,
                                 double surfaceTilt = 90,
                                 bool simulateVaporFluxEnergy = true);
    };

    ////////////////////////////////////////////////////////
    /// ASHRAEOutsideConvectionBC
    ////////////////////////////////////////////////////////

    //! \brief Convection boundary condition with ASHRAE outside convection algorithm
    class ASHRAEOutsideConvectionBC : public IConvectionBC
    {
    public:
        //! \brief Construction of ASHRAE outside convection algorithm
        //!
        //! \param index1 Node 1 index
        //! \param index2 Node 2 index
        //! \param coeff Coefficients for that are variable within timestep
        //! \param simulateVaporFluxEnergy Flag to indicate whether or not to include energy from
        //! vapor and water flux.
        ASHRAEOutsideConvectionBC(size_t index1,
                                  size_t index2,
                                  const ASHRAEOutsideCoefficients & coeff,
                                  bool simulateVaporFluxEnergy = true);
    };

    ////////////////////////////////////////////////////////
    /// YazdanianKlemsConvectionBC
    ////////////////////////////////////////////////////////

    //! \brief Convection boundary condition with Yazdanian-Klems convection algorithm.
    class YazdanianKlemsConvectionBC : public IConvectionBC
    {
    public:
        //! \brief Construction of Yazdanian-Klems outside convection algorithm
        //!
        //! \param index1 Node 1 index
        //! \param index2 Node 2 index
        //! \param coeff Coefficients for Yazdanian-Klems convection algorithm that are variable
        //! within timestep
        //! \param simulateVaporFluxEnergy Flag to indicate whether or not to
        //! include energy from vapor and water flux.
        YazdanianKlemsConvectionBC(size_t index1,
                                   size_t index2,
                                   const YazdanianKlemsCoefficients & coeff,
                                   bool simulateVaporFluxEnergy = true);
    };

    ////////////////////////////////////////////////////////
    /// KimuraConvectionBC
    ////////////////////////////////////////////////////////

    //! \brief Convection boundary condition with Kimura convection algorithm.
    class KimuraConvectionBC : public IConvectionBC
    {
    public:
        //! \brief Construction of Kimura convection algorithm
        //!
        //! \param index1 Node 1 index
        //! \param index2 Node 2 index
        //! \param coeff Coefficients Coefficients for Kimura convection algorithm that are variable
        //! within timestep
        //! \param simulateVaporFluxEnergy Flag to indicate whether or not to include energy from
        //! vapor and water flux.
        KimuraConvectionBC(size_t index1,
                           size_t index2,
                           const KimuraCoefficients & coeff,
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
        [[nodiscard]] std::vector<double> R_Vector() const override;

        //! Function that calculates matrix.
        [[nodiscard]] SquareMatrix H_Matrix() const override;

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
        [[nodiscard]] std::vector<double> R_Vector() const override;

        //! Function that calculates matrix.
        [[nodiscard]] SquareMatrix H_Matrix() const override;

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
}   // namespace HygroThermFEM
