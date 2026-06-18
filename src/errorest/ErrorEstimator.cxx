#include <algorithm>
#include <cmath>
#include <cstddef>
#include <ranges>
#include <vector>

#include <SquareMatrix.hpp>

#include "lbnl/algorithm.hxx"
#include "lbnl/expected.hxx"

#include "ErrorEstimator.hxx"
#include "FeKernels.hxx"

namespace lbnl::errorest
{
    namespace
    {
        using detail::applyTensor;
        using detail::dot2;
        using detail::dot4;
        using detail::gaussRule;
        using detail::jacobianDet;
        using detail::recoveryBasis;
        using detail::shapeFunctions;
        using FenestrationCommon::LUFactor;
        using FenestrationCommon::SquareMatrix;

        //! Number of monomials in the bilinear recovery basis [1, x, y, x*y].
        constexpr std::size_t kOrder = 4;

        bool allFinite(const std::array<double, 4> & vec)
        {
            return std::ranges::all_of(vec, [](const double val) { return std::isfinite(val); });
        }

        //! Add one Gauss point's contribution to a vertex patch's SPR normal equations.
        //! Right-hand-side columns 0 and 1 carry the two flux components.
        void accumulatePatch(SquareMatrix & matA,
                             SquareMatrix & rhsB,
                             const std::array<double, 4> & bas,
                             const Flux & sig)
        {
            for (std::size_t row = 0; row < kOrder; ++row)
            {
                for (std::size_t col = 0; col < kOrder; ++col)
                {
                    matA(row, col) += bas[row] * bas[col];
                }
                rhsB(row, 0) += bas[row] * sig[0];
                rhsB(row, 1) += bas[row] * sig[1];
            }
        }

        std::array<double, 4> column(const SquareMatrix & mat, const std::size_t col)
        {
            return {mat(0, col), mat(1, col), mat(2, col), mat(3, col)};
        }

        //! Superconvergent patch recovery: a smoothed nodal flux at every touched node.
        lbnl::ExpectedExt<std::vector<Flux>, Error> recoverNodalFlux(const Input & inp)
        {
            const std::size_t nodeCount = inp.nodes.size();
            std::vector<SquareMatrix> matA(nodeCount, SquareMatrix(kOrder));
            std::vector<SquareMatrix> rhsB(nodeCount, SquareMatrix(kOrder));
            std::vector<char> touched(nodeCount, 0);

            for (const auto & ele : inp.elements)
            {
                for (std::size_t pnt = 0; pnt < kOrder; ++pnt)
                {
                    const auto bas = recoveryBasis(ele.gaussPoints[pnt]);
                    const auto & sig = ele.flux[pnt];
                    for (int vtx = 0; vtx < ele.vertexCount; ++vtx)
                    {
                        const auto nid = static_cast<std::size_t>(ele.vertexIds[vtx]);
                        touched[nid] = 1;
                        accumulatePatch(matA[nid], rhsB[nid], bas, sig);
                    }
                }
            }

            std::vector<Flux> nodal(nodeCount, Flux{0.0, 0.0});
            for (std::size_t nid = 0; nid < nodeCount; ++nid)
            {
                if (touched[nid] == 0)
                {
                    continue;
                }
                const LUFactor factor(matA[nid]);
                const SquareMatrix sol = factor.solveRight(rhsB[nid]);
                const auto coefX = column(sol, 0);
                const auto coefY = column(sol, 1);
                if (!allFinite(coefX) || !allFinite(coefY))
                {
                    return lbnl::Unexpected(Error::SingularPatch);
                }
                const auto bas = recoveryBasis(inp.nodes[nid]);
                nodal[nid] = {dot4(bas, coefX), dot4(bas, coefY)};
            }
            return nodal;
        }

        //! Per-element energy and energy-norm error, integrated by the 2x2 Gauss rule.
        struct ElementError
        {
            double energy{0.0};     //!< Complementary strain energy of the FE flux.
            double error{0.0};      //!< Energy-norm error against the recovered flux.
            double relative{0.0};   //!< Relative error (fraction) for this element.
        };

        //! Gather the four corner coordinates and recovered nodal flux for an element.
        void gatherCorners(const Element & ele,
                           const Input & inp,
                           const std::vector<Flux> & nodalFlux,
                           std::array<double, 4> & xxx,
                           std::array<double, 4> & yyy,
                           std::array<Flux, 4> & star)
        {
            for (int vtx = 0; vtx < 4; ++vtx)
            {
                const int src = vtx < ele.vertexCount ? vtx : ele.vertexCount - 1;
                const auto nid = static_cast<std::size_t>(ele.vertexIds[src]);
                xxx[vtx] = inp.nodes[nid][0];
                yyy[vtx] = inp.nodes[nid][1];
                star[vtx] = nodalFlux[nid];
            }
        }

        ElementError elementError(const Element & ele,
                                 const Input & inp,
                                 const std::vector<Flux> & nodalFlux)
        {
            std::array<double, 4> xxx = {0.0, 0.0, 0.0, 0.0};
            std::array<double, 4> yyy = {0.0, 0.0, 0.0, 0.0};
            std::array<Flux, 4> star{};
            gatherCorners(ele, inp, nodalFlux, xxx, yyy, star);

            const auto rule = gaussRule();
            double energy = 0.0;
            double error = 0.0;
            for (std::size_t pnt = 0; pnt < kOrder; ++pnt)
            {
                const auto nnn = shapeFunctions(rule[pnt][0], rule[pnt][1]);
                Flux starGp = {0.0, 0.0};
                for (std::size_t vtx = 0; vtx < 4; ++vtx)
                {
                    starGp[0] += nnn[vtx] * star[vtx][0];
                    starGp[1] += nnn[vtx] * star[vtx][1];
                }
                const Flux & sig = ele.flux[pnt];
                const Flux err = {starGp[0] - sig[0], starGp[1] - sig[1]};
                const double dja = jacobianDet(xxx, yyy, rule[pnt][0], rule[pnt][1]);
                error += dot2(applyTensor(ele.inverseConstitutive, err), err) * dja;
                energy += dot2(applyTensor(ele.inverseConstitutive, sig), sig) * dja;
            }
            const double denom = energy + error;
            const double rel = denom == 0.0 ? 0.0 : std::sqrt(error / denom);
            return ElementError{.energy = energy, .error = error, .relative = rel};
        }
    }

    lbnl::ExpectedExt<Result, Error> estimate(const Input & inp)
    {
        if (inp.nodes.empty() || inp.elements.empty())
        {
            return lbnl::Unexpected(Error::EmptyMesh);
        }

        const auto recovered = recoverNodalFlux(inp);
        if (!recovered.has_value())
        {
            return lbnl::Unexpected(recovered.error());
        }

        Result res{};
        res.elementRelativeError.reserve(inp.elements.size());
        double totalEnergy = 0.0;
        double totalError = 0.0;
        for (const auto & ele : inp.elements)
        {
            const auto eer = elementError(ele, inp, recovered.value());
            totalEnergy += eer.energy;
            totalError += eer.error;
            res.elementRelativeError.push_back(eer.relative);
        }

        const double denom = totalEnergy + totalError;
        res.globalErrorPercent = denom == 0.0 ? 0.0 : 100.0 * std::sqrt(totalError / denom);

        const double goal = inp.targetPercent * 0.01;
        const auto indices = std::views::iota(std::size_t{0}, res.elementRelativeError.size());
        res.elementsToRefine = lbnl::transform_to_vector(
            indices | std::views::filter([&](const std::size_t idx)
                                         { return res.elementRelativeError[idx] > goal; }),
            [](const std::size_t idx) { return static_cast<int>(idx); });

        return res;
    }
}
