#pragma once

#include <memory>

#include "Elements2D.hxx"
#include "BoundaryConditions2D.hxx"
#include "FrameCavities.hxx"
#include "BoundaryConditionCoefficients.hxx"

namespace HygroThermFEM
{
    //! \brief Class to hold solution from single timestep.
    struct SingleSolution
    {
        std::vector<double> solution;   //!< Solution from single transient step
        double dTime;   //!< Timestep for which solution has been performed. Engine can adopt new
                        //!< timestep for which solution will converge.
    };

    //! \brief Interface that will keep all elements and boundary conditions together.
    //!
    //! One domain will solve single differential equation and therefore, single domain will
    //! represent thermal, moisture or pressure separately.
    class IDomain
    {
    public:
        virtual ~IDomain() = default;
        //! Domain construction. It is necessary to set up base variable that will be considered
        //! unknown.
        explicit IDomain(
          BaseVariable property,   //!< State variable which will be considered unknown.
          bool automaticUpdateOfPreviousTimestep =
            true   //!< When solver finds solution, previous timestep will be automatically updated
                   //!< by default. This should be set to false is Domain is used in outside
                   //!< iterations for solving more complex problems in which case previous timestep
                   //!< should not be updated till convergence is achieved.
        );

        //! Calculates steady state for given data
        std::vector<double> steadyState();

        //! \brief Calculates next timestep values from current (initial) values
        //! @param currentStateValues Current values of state variable or initial condition
        //! @param t_DTime Timestep in transient solution
        //! @param timestepIndex Index for current timestep. Used in variable boundary conditions
        //! case. It is defaulted to zero in case of non-variable boundary condition calculations
        //! are requested.
        SingleSolution transient(const std::vector<double> & currentStateValues,
                                 double t_DTime,
                                 size_t timestepIndex = 0);

        //! Returns flux in x and y direction
        std::vector<NodeFlux> flux() const;

        //! Adds element into domain
        virtual void createElement(
          size_t index1,                     //!< Node 1 index
          size_t index2,                     //!< Node 2 index
          size_t index3,                     //!< Node 3 index
          size_t index4,                     //!< Node 4 index
          const std::string & materialName   //!< SolidMaterial that will be assigned to the element
          ) = 0;

    protected:
        friend class MultiDomain;

        //! Forms left hand side matrix in steady state solution.
        SquareMatrix steadyStateLeftHandSide();

        //! Form right hand side vector in stead state solution.
        std::vector<double> steadyStateRightHandSide() const;

        //! \brief Forms mass, conductance and H (from boundary condition) matrices.
        //! @param t_DTime Time between two timestep for which calculates are being performed
        //! @param timestepIndex Timestep index used in case of variable timestep input boundary
        //! conditions.
        SquareMatrix transientM_K_H_Matrix(double t_DTime, size_t timestepIndex);

        //! \brief This function retrieves M*U+R vector (where U is state variable)
        //! @param t_PreviousSolution Solution from previous timestep
        //! @param t_DTime Time between two timestep for which calculates are being performed
        //! @param timestepIndex Timestep index used in case of variable timestep input boundary
        //! conditions.
        std::vector<double> transientMT_R_Vector(const std::vector<double> & t_PreviousSolution,
                                                 double t_DTime,
                                                 size_t timestepIndex);

        //! Returns if domain problem is linear.
        bool isLinear() const;

        //! Some domains require post-processing of results. Good example is
        //! moisture domain where humidity cannot go over 1.0 or lower than one.
        //! With certain set of boundary conditions and long enough time-step,
        //! solution can achieve such state and post processing should prevent it.
        virtual void postProcess(std::vector<double> & solution);

        //! \brief Calling timestep calculations
        //! @param currentStateValues Current state values from previous timestep
        //! @param t_DTime Time different for between timesteps
        //! @param timestepIndex Current timestep index used in variable boundary conditions
        std::pair<std::vector<double>, bool> transientTimestep(
          const std::vector<double> &
            currentStateValues,   //!< Current values of state variable or initial condition
          double t_DTime,         //!< Timestep in transient solution
          size_t timestepIndex);

        BaseVariable m_Property;
        ElementsLinear2D m_Elements;
        BoundaryConditions2D m_BCs;

        // Indicates if transient timestep will automatically update previous timestep solution.
        // This should be turned off if used in multidomain because previous timestep should
        // remain constant during iterations.
        bool m_AutomaticUpdatePreviousTimestep;
    };

    //! \brief Domain class for solving temperature solution.
    class ThermalDomain : public IDomain
    {
    public:
        //! Simple constructor
        ThermalDomain(bool automaticUpdatePreviousTimestep = true);

        //! \brief Creation of convection boundary condition
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param fixedBCHCCoefficients structure to hold fixed convection coefficient boundary
        //! conditions
        //! @param t_CalculateMoisture Flag on whether or not to calculate moisture
        void createConvectionBCFixedHc(size_t index1,
                                       size_t index2,
                                       const FixedBCHCCoefficients & fixedBCHCCoefficients,
                                       bool t_CalculateMoisture = true);

        //! \brief Creation of convection boundary condition
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param fixedBCHCCoefficients structure to hold fixed convection coefficients boundary
        //! conditions for transient simulations
        //! @param calculateMoisture Flag on whether or not to calculate moisture
        void createConvectionBCFixedHc(
          size_t index1,
          size_t index2,
          const std::vector<FixedBCHCCoefficients> & fixedBCHCCoefficients,
          bool calculateMoisture = true);

        //! \brief Creation of convection boundary condition
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param varHCCoeff structure to hold variable convection coefficient boundary conditions
        //! @param t_CalculateMoisture Flag on whether or not to calculate moisture
        void createConvectionBCVariableHc(size_t index1,
                                          size_t index2,
                                          const VariableBCHCCoefficients & varHCCoeff,
                                          bool t_CalculateMoisture = true);

        //! \brief Creation of convection boundary condition
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param varHCCoeff structure to hold variable convection coefficients boundary
        //! conditions for transient simulations
        //! @param t_CalculateMoisture Flag on whether or not to calculate moisture
        void createConvectionBCVariableHc(size_t index1,
                                          size_t index2,
                                          const std::vector<VariableBCHCCoefficients> & varHCCoeff,
                                          bool t_CalculateMoisture = true);

        //! \brief Creation of temperature boundary condition
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param t_Temp1 Constant temperature at node 1
        //! @param t_Temp2 Constant temperature at node 2
        void createTemperatureBC(size_t index1, size_t index2, double t_Temp1, double t_Temp2);


        //! \brief Sets fixed temperature boundary conditions
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param temp Set of node temperatures for every given timstep (each node can have
        //! different temperature).
        void createTemperatureBC(size_t index1,
                                 size_t index2,
                                 const std::vector<ConstantBCTemperatures> & temp);

        //! \brief Creation of temperature boundary condition
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param t_Temp Constant temperature in both nodes
        void createTemperatureBC(size_t index1, size_t index2, double t_Temp);

        //! \brief Creation of temperature boundary condition
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param t_Temp Constant temperature at each timestep in both nodes
        void createTemperatureBC(size_t index1, size_t index2, std::vector<double> t_Temp);

        //! \brief Creation of flux boundary condition
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param t_Flux Constant flux in both nodes
        void createFluxBC(size_t index1, size_t index2, double t_Flux);

        //! \brief Creation of black body radiation
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param t_Emissivity Emissivity at boundary condition
        //! @param t_RadiationTemperature Radiation temperature
        void createBlackBodyRadiationBC(size_t index1,
                                        size_t index2,
                                        double t_Emissivity,
                                        double t_RadiationTemperature);

        //! \brief Creation of linearized radiation boundary condition
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param linearRadBC Structure that holds coefficients for linearized radiation
        void createSimplifiedRadiationBC(size_t index1,
                                         size_t index2,
                                         const LinearizedRadiationBCCoefficients & linearRadBC);

        //! \brief Creates and adds element into domain.
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param index3 Node 3 index
        //! @param index4 Node 4 index
        //! @param materialName Material name assigned to the element
        void createElement(size_t index1,
                           size_t index2,
                           size_t index3,
                           size_t index4,
                           const std::string & materialName) override;

    protected:
        //! Storage for frame cavities recalculation
        std::unique_ptr<EquivalentFrameCavities> frameCavities;

        void postProcess(std::vector<double> & solution) override;
    };

    //! \brief Domain class for solving humidity distribution.
    class MoistureDomain : public IDomain
    {
    public:
        //! Simple constructor
        MoistureDomain(bool automaticUpdatePreviousTimestep = true);

        //! \brief Creation of convection boundary condition
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param varHCCoeff structure to hold variable convection coefficient boundary conditions
        void createMoistureBCVariableHc(size_t index1,
                                        size_t index2,
                                        const VariableBCHCCoefficients & varHCCoeff);

        //! \brief Creation of convection boundary condition
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param varCoeff structure to hold fixed convection coefficient boundary
        //! conditions for every timestep
        void createMoistureBCVariableHc(size_t index1,
                                        size_t index2,
                                        const std::vector<VariableBCHCCoefficients> & varCoeff);

        //! \brief Creation of convection boundary condition
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param fixedBchcCoefficients structure to hold fixed convection coefficient boundary
        //! conditions
        void createMoistureBCFixedHc(size_t index1,
                                     size_t index2,
                                     const FixedBCHCCoefficients & fixedBchcCoefficients);

        //! \brief Creation of convection boundary condition
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param fixedBchcCoefficients structure to hold fixed convection coefficient boundary
        //! conditions for every timestep
        void
          createMoistureBCFixedHc(size_t index1,
                                  size_t index2,
                                  const std::vector<FixedBCHCCoefficients> & fixedBchcCoefficients);

        //! \brief Creates and adds element into domain.
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param index3 Node 3 index
        //! @param index4 Node 4 index
        //! @param materialName Material name assigned to the element
        virtual void createElement(size_t index1,
                                   size_t index2,
                                   size_t index3,
                                   size_t index4,
                                   const std::string & materialName) override;

    protected:
        void postProcess(std::vector<double> & solution) override;
    };

}   // namespace HygroThermFEM
