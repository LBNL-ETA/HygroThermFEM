#pragma once

#include <memory>

#include "Elements2D.hxx"
#include "BoundaryConditions2D.hxx"
#include "GasCavities.hxx"
#include "BoundaryConditionCoefficients.hxx"
#include "TimestepNotifier.hxx"
#include "TimestepObserver.hxx"
#include "Solver.hxx"

namespace HygroThermFEM
{
    //! \brief Interface that will keep all elements and boundary conditions together.
    //!
    //! One domain will solve single differential equation and therefore, single domain will
    //! represent thermal, moisture or pressure separately.
    struct SingleDomain : public Timesteps::TimestepNotifier
    {
        virtual ~SingleDomain() = default;
        //! Domain construction. It is necessary to set up base variable that will be considered
        //! unknown.
        //! @param property State variable which will be considered unknown.
        //! @param automaticUpdateOfPreviousTimestep When solver finds solution, previous timestep
        //! will be automatically updated by default. This should be set to false is Domain is used
        //! in outside iterations for solving more complex problems in which case previous timestep
        //! should not be updated till convergence is achieved.
        explicit SingleDomain(
          BaseVariable property,
          bool automaticUpdateOfPreviousTimestep = true);

        //! Returns flux in x and y direction
        [[nodiscard]] std::vector<NodeFlux> flux() const;

        //! Adds element into domain
        virtual void createElement(
          size_t index1,                     //!< Node 1 index
          size_t index2,                     //!< Node 2 index
          size_t index3,                     //!< Node 3 index
          size_t index4,                     //!< Node 4 index
          const std::string & materialName   //!< SolidMaterial that will be assigned to the element
          ) = 0;

        //! \brief Sets new gravity vector and performs new calculations
        //!
        //! @param gravityVector Direction of gravity
        void setGravityVector(const FenestrationCommon::GravityVector & gravityVector);

        //! Some domains require post-processing of results. Good example is
        //! moisture domain where humidity cannot go over 1.0 or lower than one.
        //! With certain set of boundary conditions and long enough time-step,
        //! solution can achieve such state and post processing should prevent it.
        virtual void postProcess(std::vector<double> & solution);

        BaseVariable m_Property;
        ElementsLinear2D m_Elements;
        BoundaryConditions2D m_BCs;

        FenestrationCommon::GravityVector m_GravityVector{0, -1, 0};

        //! Storage for gas cavities recalculation
        std::unique_ptr<EquivalentGasCavities> gasCavities;

        // Indicates if transient timestep will automatically update previous timestep solution.
        // This should be turned off if used in multidomain because previous timestep should
        // remain constant during iterations.
        bool m_AutomaticUpdatePreviousTimestep;
    };    

}   // namespace HygroThermFEM
