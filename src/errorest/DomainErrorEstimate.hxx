#pragma once

#include "ErrorEstimator.hxx"

namespace HygroThermFEM
{
    class IDomain;
}

//! Bridge between a solved HygroThermFEM domain and the field-agnostic recovery-based error
//! estimator in lbnl::errorest. The estimator itself never references the domain; this adapter
//! is the one place that reads a solved domain's elements and assembles the estimator Input.
namespace HygroThermFEM::errorest
{
    //! The estimator's plain-data Input (per-element geometry + FE flux) together with its Result.
    //! The Input is returned so the caller can serialize the full per-element data -- e.g. the
    //! legacy mesher's ".eedump" -- without rebuilding it from the domain.
    struct DomainEstimate
    {
        lbnl::errorest::Input input{};
        lbnl::errorest::Result result{};
    };

    //! \brief Run the recovery-based (Zienkiewicz-Zhu) error estimator on a solved domain.
    //!
    //! Assembles lbnl::errorest::Input from the domain's elements: the global Gauss-point
    //! coordinates, the FE flux (sigma = -k * grad) at those points, and the per-element inverse
    //! conductivity that defines the energy norm. Elements are grouped for recovery by material so
    //! patches do not smooth flux across a material interface.
    //!
    //! The domain must already be solved (nodal state populated). The estimator can only fail on a
    //! degenerate or empty mesh; that is reported as an all-zero Result so the caller treats it as
    //! "no actionable error estimate" and stops refining rather than looping forever.
    //!
    //! \param domain Solved domain whose discretization error is estimated.
    //! \param targetPercent Accuracy target [percent]; elements above it appear in
    //!        Result::elementsToRefine.
    //! \return The Input that was fed to the estimator and the Result it produced.
    [[nodiscard]] DomainEstimate estimateError(const IDomain & domain, double targetPercent);
}
