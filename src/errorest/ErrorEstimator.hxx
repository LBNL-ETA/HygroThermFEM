#pragma once

#include <array>
#include <vector>

#include "lbnl/expected.hxx"

//! Recovery-based (Zienkiewicz-Zhu) a-posteriori error estimator.
//!
//! Field-agnostic by design: the estimator never references temperature,
//! moisture, or any specific physics. The only physical input is the per-element
//! inverse constitutive tensor that defines the energy norm. See the methodology
//! document `D:\Documents\Error Estimator\error_estimation_methodology.typ`.
namespace lbnl::errorest
{
    using Point = std::array<double, 2>;    //!< Global (x, y) coordinate.
    using Flux = std::array<double, 2>;     //!< Flux vector (sigma_x, sigma_y).
    using Tensor = std::array<double, 4>;   //!< Row-major 2x2 [d00, d01, d10, d11].

    //! One finite element's contribution to the estimate, as plain data.
    struct Element
    {
        std::array<Point, 4> gaussPoints{}; //!< Global coords of the 2x2 Gauss points.
        std::array<Flux, 4> flux{};         //!< FE flux sigma_h at each Gauss point.
        Tensor inverseConstitutive{};       //!< D^-1 used to weight the energy norm.
        std::array<int, 4> vertexIds{};     //!< Node indices; a triangle duplicates the last.
        int vertexCount{4};                 //!< 3 (triangle) or 4 (quad).
    };

    //! Complete, file-free input to the estimator.
    struct Input
    {
        std::vector<Point> nodes{};         //!< Node coordinates, indexed by Element::vertexIds.
        std::vector<Element> elements{};    //!< One entry per finite element.
        double targetPercent{0.0};          //!< Accuracy target, percent.
    };

    //! Plain-data result; the caller persists or acts on it.
    struct Result
    {
        double globalErrorPercent{0.0};         //!< Global relative error, percent.
        std::vector<double> elementRelativeError{};  //!< Per-element relative error, fraction.
        std::vector<int> elementsToRefine{};    //!< Indices of elements above the target.
    };

    enum class Error
    {
        EmptyMesh,      //!< No nodes or no elements supplied.
        SingularPatch,  //!< Degenerate geometry; an SPR patch was rank deficient.
    };

    //! Estimate the discretization error of an FE solution.
    [[nodiscard]] lbnl::ExpectedExt<Result, Error> estimate(const Input & inp);
}
