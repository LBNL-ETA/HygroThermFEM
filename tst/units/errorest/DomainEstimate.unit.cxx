#include <array>

#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "DomainErrorEstimate.hxx"
#include "TestMaterials.hxx"

//! Adapter-level tests for HygroThermFEM::errorest::estimateError. They build a domain whose nodal
//! field is prescribed to a known function and then run the recovery-based estimator on it, so they
//! cover the whole bridge deterministically: per-element Gauss-point flux, the global Gauss
//! coordinates, and the inverse conductivity that defines the energy norm. The end-to-end solve is
//! exercised by the THERM steady-state path; the estimator math is covered by Estimator.unit /
//! Regression. Driving the field directly keeps these tests independent of solver/BC behaviour.

namespace
{
    //! Builds a 2x2 grid of bilinear quads on [0,2] x [0,2] and sets every node's temperature to
    //! the supplied sampled field (node number n carries field[n - 1]). The domain is populated in
    //! place because MultiDomain is neither copyable nor movable.
    void buildGrid(HygroThermFEM::MultiDomain & multiDomain,
                   const std::array<double, 9> & field)
    {
        auto & nodes = multiDomain.nodes();
        const auto node = [&](size_t index, double xpos, double ypos) {
            nodes.createNode({.index = index,
                              .x = xpos,
                              .y = ypos,
                              .state = HygroThermFEM::State{.temperature = field[index - 1]}});
        };
        node(1, 0, 0);
        node(2, 1, 0);
        node(3, 2, 0);
        node(4, 0, 1);
        node(5, 1, 1);
        node(6, 2, 1);
        node(7, 0, 2);
        node(8, 1, 2);
        node(9, 2, 2);

        const auto & material = multiDomain.materials().createSolidMaterial(TestHelper::TestMaterial());
        const auto & mat = material.name();

        // Counter-clockwise winding: bottom-left, bottom-right, top-right, top-left.
        multiDomain.createElement({.node1 = 1, .node2 = 2, .node3 = 5, .node4 = 4, .material = mat});
        multiDomain.createElement({.node1 = 2, .node2 = 3, .node3 = 6, .node4 = 5, .material = mat});
        multiDomain.createElement({.node1 = 4, .node2 = 5, .node3 = 8, .node4 = 7, .material = mat});
        multiDomain.createElement({.node1 = 5, .node2 = 6, .node3 = 9, .node4 = 8, .material = mat});
    }
}

TEST(DomainEstimate, LinearFieldHasNearZeroError)
{
    SCOPED_TRACE("A linear field has constant gradient, so the FE flux is constant and the "
                 "recovered flux matches it exactly -- the estimated error must be ~0.");

    // T(x, y) = 10 + 5x + 3y sampled at the nine nodes.
    const std::array<double, 9> linear{10, 15, 20, 13, 18, 23, 16, 21, 26};
    HygroThermFEM::MultiDomain multiDomain({.performThermal = true, .performMoisture = false});
    buildGrid(multiDomain, linear);

    const auto result = HygroThermFEM::errorest::estimateError(multiDomain.thermal(), 2.0);

    EXPECT_LT(result.globalErrorPercent, 1e-6);
    EXPECT_TRUE(result.elementsToRefine.empty());
}

TEST(DomainEstimate, NonLinearFieldHasErrorAndFlagsElements)
{
    SCOPED_TRACE("T = x^2 + y^2 has a varying gradient the coarse mesh cannot reproduce, so the "
                 "estimator must report a positive error and flag elements for refinement.");

    // T(x, y) = x^2 + y^2 sampled at the nine nodes.
    const std::array<double, 9> quadratic{0, 1, 4, 1, 2, 5, 4, 5, 8};
    HygroThermFEM::MultiDomain multiDomain({.performThermal = true, .performMoisture = false});
    buildGrid(multiDomain, quadratic);

    const auto result = HygroThermFEM::errorest::estimateError(multiDomain.thermal(), 0.0);

    EXPECT_GT(result.globalErrorPercent, 0.0);
    EXPECT_FALSE(result.elementsToRefine.empty());
}
