#pragma once

#include <memory>
#include <optional>

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

    //! \brief Interface that will keep all elements and boundary conditions together.
    //!
    //! One domain will solve single differential equation and therefore, single domain will
    //! represent thermal, moisture or pressure separately.
    struct SingleDomain : public Timesteps::TimestepNotifier
    {
        //! Domain construction. It is necessary to set up base variable that will be considered
        //! unknown.
        //! @param type Type of domain that will be solved for (Thermal or Moisture).
        explicit SingleDomain(DomainType type);

        //! Returns flux in x and y direction
        [[nodiscard]] std::vector<NodeFlux> flux() const;

        //! \brief Sets new gravity vector and performs new calculations
        //!
        //! @param gravityVector Direction of gravity
        void setGravityVector(const FenestrationCommon::GravityVector & gravityVector);

        DomainType domainType;

        ElementsLinear2D elements;
        BoundaryConditions2D boundaryConditions;

        FenestrationCommon::GravityVector gravityVector{0, -1, 0};

        //! Storage for gas cavities recalculation
        std::optional<EquivalentGasCavities> gasCavities{std::nullopt};
    };

    BaseVariable baseVariable(SingleDomain & domain);

}   // namespace HygroThermFEM
