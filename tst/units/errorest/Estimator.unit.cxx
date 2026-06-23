#include <array>

#include <gtest/gtest.h>

#include "ErrorEstimator.hxx"
#include "FeKernels.hxx"

using namespace lbnl::errorest;

namespace
{
    //! A single unit-square quad whose Gauss-point flux is given by `sampler`.
    //! Nodes are the corners (0,0),(1,0),(1,1),(0,1); the Gauss points are mapped
    //! from the reference square via x = (r+1)/2, y = (s+1)/2.
    template<typename Sampler>
    Input singleSquare(Sampler sampler)
    {
        const auto rule = detail::gaussRule();
        Element ele{};
        for (std::size_t pnt = 0; pnt < 4; ++pnt)
        {
            const Point pos = {(rule[pnt][0] + 1.0) / 2.0, (rule[pnt][1] + 1.0) / 2.0};
            ele.gaussPoints[pnt] = pos;
            ele.flux[pnt] = sampler(pos);
        }
        ele.inverseConstitutive = {1.0, 0.0, 0.0, 1.0};
        ele.vertexIds = {0, 1, 2, 3};   // distinct 4th node => quad (isTriangle == false)

        Input inp{};
        inp.nodes = {{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}};
        inp.elements = {ele};
        inp.targetPercent = 2.0;
        return inp;
    }
}

TEST(Estimator, EmptyMeshIsError)
{
    const auto res = estimate(Input{});
    ASSERT_FALSE(res.has_value());
    EXPECT_EQ(res.error(), Error::EmptyMesh);
}

TEST(Estimator, MixedSubdomainsIsError)
{
    // Region topology is all-or-nothing: one element carrying a subdomain while
    // another does not is malformed input and must be reported, not papered over.
    auto inp = singleSquare([](const Point &) { return Flux{1.0, 0.0}; });
    inp.elements[0].subdomain = 0;          // assigned
    Element unassigned = inp.elements[0];
    unassigned.subdomain.reset();           // the malformed case
    inp.elements.push_back(unassigned);
    const auto res = estimate(inp);
    ASSERT_FALSE(res.has_value());
    EXPECT_EQ(res.error(), Error::UnassignedSubdomain);
}

TEST(Estimator, ConstantFluxHasZeroError)
{
    // A constant flux is recovered exactly, so the estimated error must vanish.
    const auto inp = singleSquare([](const Point &) { return Flux{5.0, -3.0}; });
    const auto res = estimate(inp);
    ASSERT_TRUE(res.has_value());
    EXPECT_NEAR(res->globalErrorPercent, 0.0, 1e-9);
    EXPECT_TRUE(res->elementsToRefine.empty());
}

TEST(Estimator, LinearFluxInSpanHasZeroError)
{
    // A flux field inside the bilinear recovery span [1, x, y, x*y] is also
    // recovered exactly on a single element, so the error is still zero.
    const auto sampler = [](const Point & pos)
    {
        return Flux{1.0 + 2.0 * pos[0] - 0.5 * pos[1], 3.0 * pos[0] * pos[1]};
    };
    const auto res = estimate(singleSquare(sampler));
    ASSERT_TRUE(res.has_value());
    EXPECT_NEAR(res->globalErrorPercent, 0.0, 1e-7);
}
