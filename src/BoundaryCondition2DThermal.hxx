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
        //! @param nodePool Reference to NodePool for node lookup
        //! @param index1 node1 index
        //! @param index2 node2 index
        //! @param fixedBCHCCoefficients Structure that holds necessary coefficients in other to
        //! create constant film coefficient boundary condition
        //! @param simulateVaporFluxEnergy Need to be true if we want to include vapor flux energy
        ConstantConvectionBC(Nodes & nodePool,
                             size_t index1,
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
        //! @param nodePool Reference to NodePool for node lookup
        //! @param index1 node1 index
        //! @param index2 node2 index
        //! @param varHCCoeff Necessary coefficients that will be used to create TARP heat flow
        //! convection algorithm
        //! @param surfaceTilt Surface tilt at the boundary
        //! @param simulateVaporFluxEnergy Indicates whether or not moisture will be calculated at
        //! the boundary
        ThermalTARPConvectionBC(Nodes & nodePool,
                                size_t index1,
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
        //! @param nodePool Reference to NodePool for node lookup
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param coeff Timestep coefficients that are used in ASHRAE Inside convection algorithm
        //! @param surfaceHeight Surface height for which convection coefficient has been calculated
        //! [meters]
        //! @param surfaceTilt Surface tilt at the boundary [degrees]
        //! @param simulateVaporFluxEnergy Indicates whether energy from vapor flow should be
        //! included in calculations.
        ASHRAEInsideConvectionBC(Nodes & nodePool,
                                 size_t index1,
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
        //! \param nodePool Reference to NodePool for node lookup
        //! \param index1 Node 1 index
        //! \param index2 Node 2 index
        //! \param coeff Coefficients for that are variable within timestep
        //! \param simulateVaporFluxEnergy Flag to indicate whether or not to include energy from
        //! vapor and water flux.
        ASHRAEOutsideConvectionBC(Nodes & nodePool,
                                  size_t index1,
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
        //! \param nodePool Reference to NodePool for node lookup
        //! \param index1 Node 1 index
        //! \param index2 Node 2 index
        //! \param coeff Coefficients for Yazdanian-Klems convection algorithm that are variable
        //! within timestep
        //! \param simulateVaporFluxEnergy Flag to indicate whether or not to
        //! include energy from vapor and water flux.
        YazdanianKlemsConvectionBC(Nodes & nodePool,
                                   size_t index1,
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
        //! \param nodePool Reference to NodePool for node lookup
        //! \param index1 Node 1 index
        //! \param index2 Node 2 index
        //! \param coeff Coefficients Coefficients for Kimura convection algorithm that are variable
        //! within timestep
        //! \param simulateVaporFluxEnergy Flag to indicate whether or not to include energy from
        //! vapor and water flux.
        KimuraConvectionBC(Nodes & nodePool,
                           size_t index1,
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
        TemperatureBC(Nodes & nodePool,        //!< Reference to NodePool for node lookup
                      size_t index1,              //!< Node 1 index
                      size_t index2,              //!< Node 2 index
                      double t_NodeTemperatures   //!< Temperature at both nodes
        );

        //! Construction of temperature boundary condition with two nodes and temperature at each of
        //! them
        TemperatureBC(Nodes & nodePool,  //!< Reference to NodePool for node lookup
                      size_t index1,        //!< Node 1 index
                      size_t index2,        //!< Node 2 index
                      double t_Temp1,       //!< Temperature at node 1
                      double t_Temp2        //!< Temperature at node 2
        );
    private:
        inline static const double hugeFilmCoefficient{1e18};
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
        FluxBC(Nodes & nodePool,  //!< Reference to NodePool for node lookup
               size_t index1,        //!< Node 1 index
               size_t index2,        //!< Node 2 index
               double t_Flux         //!< Constant flux at both nodes
        );

        //! Function that calculates right hand side vector.
        [[nodiscard]] std::vector<double> R_Vector() const override;

        //! Function that calculates matrix.
        [[nodiscard]] SquareMatrix H_Matrix() const override;

    private:
        double m_Flux;
    };

    ///////////////////////////////////////////////////////
    /// IRadiationBC - Base class for radiation BCs
    ///////////////////////////////////////////////////////

    //! \brief Base class for radiation boundary conditions.
    //!
    //! Provides common implementation of R_Vector() and H_Matrix() for radiation-based
    //! boundary conditions. Derived classes implement radiationCoefficients() to provide
    //! the specific radiation calculation method.
    class IRadiationBC : public IBCLinear2D
    {
    public:
        virtual ~IRadiationBC() = default;

        //! Function that calculates right hand side vector.
        [[nodiscard]] std::vector<double> R_Vector() const override;

        //! Function that calculates matrix.
        [[nodiscard]] SquareMatrix H_Matrix() const override;

    protected:
        //! Construction of radiation boundary condition base.
        IRadiationBC(Nodes & nodePool,            //!< Reference to NodePool for node lookup
                     size_t index1,                  //!< Node 1 index
                     size_t index2,                  //!< Node 2 index
                     double radiationTemperature,    //!< Outside radiation temperature in celsius.
                     bool isLinear = true            //!< Whether the BC is linear
        );

        //! Returns radiation coefficients for each boundary node.
        //! Derived classes implement specific calculation methods.
        [[nodiscard]] virtual std::vector<double> radiationCoefficients() const = 0;

        double m_RadiationTemperature;
    };

    ///////////////////////////////////////////////////////
    /// BlackBodyRadiationBC
    ///////////////////////////////////////////////////////

    //! \brief Simple black body radiation boundary condition.
    //!
    //! Boundary condition is used only in thermal module. Uses Stefan-Boltzmann law
    //! to calculate radiative heat transfer coefficient based on current temperatures.
    class BlackBodyRadiationBC : public IRadiationBC
    {
    public:
        //! Black body radiation boundary condition with two nodes, emissivity and radiation
        //! temperature
        BlackBodyRadiationBC(
          Nodes & nodePool,            //!< Reference to NodePool for node lookup
          size_t index1,                  //!< Node 1 index
          size_t index2,                  //!< Node 2 index
          double emissivity,              //!< Boundary condition surface emissivity at both nodes.
          double radiationTemperature     //!< Outside radiation temperature in celsius.
        );

    protected:
        //! Radiative coefficient calculated using Stefan-Boltzmann law based on current temperatures
        [[nodiscard]] std::vector<double> radiationCoefficients() const override;

    private:
        double m_Emissivity;
    };

    ///////////////////////////////////////////////////////
    /// LinearizedRadiationBC
    ///////////////////////////////////////////////////////

    //! \brief Linearized radiation boundary condition.
    //!
    //! Uses a pre-calculated constant radiation coefficient instead of
    //! computing it from Stefan-Boltzmann law at each iteration.
    class LinearizedRadiationBC : public IRadiationBC
    {
    public:
        LinearizedRadiationBC(Nodes & nodePool,
                              size_t index1,
                              size_t index2,
                              const LinearizedRadiationBCCoefficients & linearRadBC);

    protected:
        //! Returns constant radiation coefficient for each boundary node
        [[nodiscard]] std::vector<double> radiationCoefficients() const override;

    private:
        double m_RadiationCoefficient;
    };
}   // namespace HygroThermFEM
