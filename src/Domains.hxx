#pragma once

#include <memory>
#include "Elements2D.hxx"
#include "BoundaryConditions2D.hxx"
#include "Elements2D.hxx"

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
        //! Domain construction. It is necessary to set up base variable that will be considered
        //! unknown.
        explicit IDomain(
          const BaseVariable property,   //!< State variable which will be considered unknown.
          bool automaticUpdateOfPreviousTimestep =
            true   //!< When solver finds solution, previous timestep will be automatically updated
                   //!< by default. This should be set to false is Domain is used in outside
                   //!< iterations for solving more complex problems in which case previous timestep
                   //!< should not be updated till convergence is achieved.
        );

        //! Calculates steady state for given data
        std::vector<double> steadyState();

        //! Calculates next timestep values from current (initial) values
        SingleSolution transient(
          const std::vector<double> &
            currentStateValues,   //!< Current values of state variable or initial condition
          double t_DTime          //!< Timestep in transient solution
        );

        //! Returns flux in x and y direction
        std::vector<NodeFlux> flux() const;

        //! Adds element into domain
        virtual void createElement(
          size_t index1,                     //!< Node 1 index
          size_t index2,                     //!< Node 2 index
          size_t index3,                     //!< Node 3 index
          size_t index4,                     //!< Node 4 index
          const std::string & materialName   //!< Material that will be assigned to the element
          ) = 0;

    protected:
        friend class MultiDomain;

        //! Helper class that will update property at nodes from new timestep values.
        void updateNodeValues(const std::vector<double> & values,
                              const BaseVariable property,
                              bool updatePreviousValues = true);

        //! Forms left hand side matrix in steady state solution.
        FenestrationCommon::SquareMatrix steadyStateLeftHandSide();

        //! Form right hand side vector in stead state solution.
        std::vector<double> steadyStateRightHandSide() const;

        //! Forms mass, conductance and H (from boundary condition) matrices.
        FenestrationCommon::SquareMatrix transientM_K_H_Matrix(const double t_DTime);

        //! This function retrieves M*U+R vector (where U is state variable)
        std::vector<double> transientMT_R_Vector(const std::vector<double> & t_PreviousSolution,
                                                 const double t_DTime);

        //! Returns if domain problem is linear.
        bool isLinear() const;

        //! Some domains require post-processing of results. Good example is
        //! moisture domain where humidity cannot go over 1.0 or lower than one.
        //! With certain set of boundary conditions and long enough time-step,
        //! solution can achieve such state and post processing should prevent it.
        virtual void postProcess(std::vector<double> & solution) const;

        std::pair<std::vector<double>, bool> transientTimestep(
          const std::vector<double> &
            currentStateValues,   //!< Current values of state variable or initial condition
          double t_DTime          //!< Timestep in transient solution
        );

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

        //! Creation of convection boundary condition
        void createConvectionBCFixedHc(size_t index1, size_t index2, double t_AirTemperature,
                                       double t_ConvectionCoefficient);

        //! Creation of convection boundary condition
        void createConvectionBCVariableHc(size_t index1,                    //!< Node 1 index
                                       size_t index2,                    //!< Node 2 index
                                       double t_AirTemperature           //!< Outside air temperature
        );

        //! Creation of temperature boundary condition
        void createTemperatureBC(size_t index1,    //!< Node 1 index
                                 size_t index2,    //!< Node 2 index
                                 double t_Temp1,   //!< Constant temperature at node 1
                                 double t_Temp2    //!< Constant temperature at node 2
        );

        //! Creation of temperature boundary condition
        void createTemperatureBC(size_t index1,   //!< Node 1 index
                                 size_t index2,   //!< Node 2 index
                                 double t_Temp    //!< Constant temperature in both nodes.
        );

        //! Creation of flux boundary condition
        void createFluxBC(size_t index1,   //!< Node 1 index
                          size_t index2,   //!< Node 2 index
                          double t_Flux    //!< Constant flux in both nodes.
        );

        //! Creation of black body radiation
        void createBlackBodyRadiationBC(size_t index1,         //!< Node 1 index
                                        size_t index2,         //!< Node 2 index
                                        double t_Emissivity,   //! Emissivity at boundary condition
                                        double t_RadiationTemperature   //! Radiation temperature
        );

        //! Creates and adds element into domain.
        virtual void createElement(size_t index1,                     //!< Node 1 index
                                   size_t index2,                     //!< Node 2 index
                                   size_t index3,                     //!< Node 3 index
                                   size_t index4,                     //!< Node 4 index
                                   const std::string & materialName   //!< Material name
                                   ) override;
    };

    //! \brief Domain class for solving humidity distribution.
    class MoistureDomain : public IDomain
    {
    public:
        //! Simple constructor
        MoistureDomain(bool automaticUpdatePreviousTimestep = true);

        //! Creates moisture boundary condition that will calculate convection film coefficient.
        void createMoistureBCVariableHc(const size_t index1, const size_t index2,
                                        const double t_AirHumidity, const double t_AirTemperature);

        //! Creates moisture boundary condition with fixed film coefficient
        void createMoistureBCFixedHc(const size_t index1, const size_t index2,
                                     const double t_AirTemperature,
                                     const double t_ConvectiveFilmCoefficient,
                                     const double t_AirHumidity);

        //! Creates and adds element into domain.
        virtual void createElement(size_t index1,                     //!< Node 1 index
                                   size_t index2,                     //!< Node 2 index
                                   size_t index3,                     //!< Node 3 index
                                   size_t index4,                     //!< Node 4 index
                                   const std::string & materialName   //!< Material name
                                   ) override;

    protected:
        void postProcess(std::vector<double> & solution) const override;
    };

}   // namespace HygroThermFEM
