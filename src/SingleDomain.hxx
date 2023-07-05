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
    enum class DomainType
    {
        Thermal,
        Moisture
    };

    using ElementFactory =
      std::function<std::unique_ptr<IElementLinear2D>(const size_t index1,
                                                      const size_t index2,
                                                      const size_t index3,
                                                      const size_t index4,
                                                      const std::string & materialName)>;

    using PostProcessFunc = std::function<void(std::vector<double> &)>;

    //! \brief Interface that will keep all elements and boundary conditions together.
    //!
    //! One domain will solve single differential equation and therefore, single domain will
    //! represent thermal, moisture or pressure separately.
    struct SingleDomain : public Timesteps::TimestepNotifier
    {
        //! Domain construction. It is necessary to set up base variable that will be considered
        //! unknown.
        //! @param type Type of domain that will be solved for (Thermal or Moisture).
        //! @param automaticUpdateOfPreviousTimestep When solver finds solution, previous timestep
        //! will be automatically updated by default. This should be set to false is Domain is used
        //! in outside iterations for solving more complex problems in which case previous timestep
        //! should not be updated till convergence is achieved.
        explicit SingleDomain(DomainType type);

        //! Returns flux in x and y direction
        [[nodiscard]] std::vector<NodeFlux> flux() const;

        //! Adds element into domain
        //!< Node 1 index
        //!< Node 2 index
        //!< Node 3 index
        //!< Node 4 index
        //!< SolidMaterial that will be assigned to the element
        void createElement(size_t index1,
                           size_t index2,
                           size_t index3,
                           size_t index4,
                           const std::string & materialName);

        //! \brief Sets new gravity vector and performs new calculations
        //!
        //! @param gravityVector Direction of gravity
        void setGravityVector(const FenestrationCommon::GravityVector & gravityVector);

        //! Some domains require post-processing of results. Good example is
        //! moisture domain where humidity cannot go over 1.0 or lower than one.
        //! With certain set of boundary conditions and long enough time-step,
        //! solution can achieve such state and post processing should prevent it.
        void postProcess(std::vector<double> & solution);

        DomainType domainType;

        ElementsLinear2D m_Elements;
        BoundaryConditions2D m_BCs;

        FenestrationCommon::GravityVector m_GravityVector{0, -1, 0};

        //! Storage for gas cavities recalculation
        std::unique_ptr<EquivalentGasCavities> gasCavities;
    };

    BaseVariable baseVariable(SingleDomain & domain);

}   // namespace HygroThermFEM
