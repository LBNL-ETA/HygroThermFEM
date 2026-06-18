#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <ranges>
#include <utility>
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

        //! A recovery patch keyed by (node, material): the SPR normal equations for one
        //! material's contribution at one node. Keying by material keeps patches from
        //! smoothing flux across a material interface (k jump), mirroring the legacy
        //! esterror, which recovers per vertex-use (calcSmoothed.c: VU_faces).
        using PatchKey = std::pair<std::size_t, int>;

        struct Patch
        {
            SquareMatrix matA;
            SquareMatrix rhsB;

            Patch() : matA(kOrder), rhsB(kOrder) {}
        };

        bool allFinite(const std::array<double, 4> & vec)
        {
            return std::ranges::all_of(vec, [](const double val) { return std::isfinite(val); });
        }

        //! Add one Gauss point's contribution to a patch's SPR normal equations.
        //! Right-hand-side columns 0 and 1 carry the two flux components.
        void accumulatePatch(Patch & patch,
                             const std::array<double, 4> & bas,
                             const Flux & sig)
        {
            for (std::size_t row = 0; row < kOrder; ++row)
            {
                for (std::size_t col = 0; col < kOrder; ++col)
                {
                    patch.matA(row, col) += bas[row] * bas[col];
                }
                patch.rhsB(row, 0) += bas[row] * sig[0];
                patch.rhsB(row, 1) += bas[row] * sig[1];
            }
        }

        std::array<double, 4> column(const SquareMatrix & mat, const std::size_t col)
        {
            return {mat(0, col), mat(1, col), mat(2, col), mat(3, col)};
        }

        //! Assign each element a small material id, one per distinct inverseConstitutive
        //! tensor, so the recovery can segregate patches by material.
        std::vector<int> materialIds(const Input & inp)
        {
            std::map<std::array<double, 4>, int> ids;
            std::vector<int> out(inp.elements.size());
            for (std::size_t idx = 0; idx < inp.elements.size(); ++idx)
            {
                const auto [pos, inserted] =
                  ids.try_emplace(inp.elements[idx].inverseConstitutive, static_cast<int>(ids.size()));
                out[idx] = pos->second;
            }
            return out;
        }

        //! Superconvergent patch recovery: a smoothed nodal flux per (node, material).
        lbnl::ExpectedExt<std::map<PatchKey, Flux>, Error> recoverNodalFlux(
          const Input & inp, const std::vector<int> & elemMat)
        {
            std::map<PatchKey, Patch> patches;
            for (std::size_t ele = 0; ele < inp.elements.size(); ++ele)
            {
                const auto & elem = inp.elements[ele];
                const int mat = elemMat[ele];
                for (std::size_t pnt = 0; pnt < kOrder; ++pnt)
                {
                    const auto bas = recoveryBasis(elem.gaussPoints[pnt]);
                    const auto & sig = elem.flux[pnt];
                    for (int vtx = 0; vtx < elem.vertexCount; ++vtx)
                    {
                        const auto nid = static_cast<std::size_t>(elem.vertexIds[vtx]);
                        accumulatePatch(patches[{nid, mat}], bas, sig);
                    }
                }
            }

            std::map<PatchKey, Flux> nodal;
            for (auto & [key, patch] : patches)
            {
                const LUFactor factor(patch.matA);
                const SquareMatrix sol = factor.solveRight(patch.rhsB);
                const auto coefX = column(sol, 0);
                const auto coefY = column(sol, 1);
                if (!allFinite(coefX) || !allFinite(coefY))
                {
                    return lbnl::Unexpected(Error::SingularPatch);
                }
                const auto bas = recoveryBasis(inp.nodes[key.first]);
                nodal[key] = {dot4(bas, coefX), dot4(bas, coefY)};
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

        //! Gather the four corner coordinates and recovered nodal flux for an element,
        //! taking the recovered flux from this element's own material patch.
        void gatherCorners(const Element & ele,
                           const int mat,
                           const Input & inp,
                           const std::map<PatchKey, Flux> & nodalFlux,
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
                star[vtx] = nodalFlux.at({nid, mat});
            }
        }

        ElementError elementError(const Element & ele,
                                 const int mat,
                                 const Input & inp,
                                 const std::map<PatchKey, Flux> & nodalFlux)
        {
            std::array<double, 4> xxx = {0.0, 0.0, 0.0, 0.0};
            std::array<double, 4> yyy = {0.0, 0.0, 0.0, 0.0};
            std::array<Flux, 4> star{};
            gatherCorners(ele, mat, inp, nodalFlux, xxx, yyy, star);

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

        const auto elemMat = materialIds(inp);
        const auto recovered = recoverNodalFlux(inp, elemMat);
        if (!recovered.has_value())
        {
            return lbnl::Unexpected(recovered.error());
        }

        Result res{};
        res.elementRelativeError.reserve(inp.elements.size());
        double totalEnergy = 0.0;
        double totalError = 0.0;
        for (std::size_t ele = 0; ele < inp.elements.size(); ++ele)
        {
            const auto eer = elementError(inp.elements[ele], elemMat[ele], inp, recovered.value());
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
