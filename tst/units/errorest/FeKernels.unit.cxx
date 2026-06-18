#include <array>
#include <cmath>

#include <gtest/gtest.h>

#include "FeKernels.hxx"

using namespace lbnl::errorest;
using namespace lbnl::errorest::detail;

TEST(FeKernels, ShapeFunctionsAreNodalAtCorners)
{
    // N_1 = 1 at the (-1,-1) corner, all others zero.
    const auto nnn = shapeFunctions(-1.0, -1.0);
    EXPECT_NEAR(nnn[0], 1.0, 1e-12);
    EXPECT_NEAR(nnn[1], 0.0, 1e-12);
    EXPECT_NEAR(nnn[2], 0.0, 1e-12);
    EXPECT_NEAR(nnn[3], 0.0, 1e-12);
}

TEST(FeKernels, ShapeFunctionsArePartitionOfUnity)
{
    const auto nnn = shapeFunctions(0.31, -0.42);
    EXPECT_NEAR(nnn[0] + nnn[1] + nnn[2] + nnn[3], 1.0, 1e-12);
}

TEST(FeKernels, RecoveryBasisIsBilinearMonomials)
{
    const auto bas = recoveryBasis(Point{2.0, 3.0});
    EXPECT_NEAR(bas[0], 1.0, 1e-12);
    EXPECT_NEAR(bas[1], 2.0, 1e-12);
    EXPECT_NEAR(bas[2], 3.0, 1e-12);
    EXPECT_NEAR(bas[3], 6.0, 1e-12);
}

TEST(FeKernels, JacobianDetOfUnitSquareIsQuarter)
{
    // Reference square has area 4; a unit physical square has area 1, so the
    // Jacobian determinant is constant 1/4 everywhere.
    const std::array<double, 4> xxx = {0.0, 1.0, 1.0, 0.0};
    const std::array<double, 4> yyy = {0.0, 0.0, 1.0, 1.0};
    EXPECT_NEAR(jacobianDet(xxx, yyy, 0.0, 0.0), 0.25, 1e-12);
    EXPECT_NEAR(jacobianDet(xxx, yyy, 0.5, -0.3), 0.25, 1e-12);
}

TEST(FeKernels, JacobianDetScalesWithArea)
{
    // A 2x3 rectangle has area 6, so detJ = 6/4 = 1.5.
    const std::array<double, 4> xxx = {0.0, 2.0, 2.0, 0.0};
    const std::array<double, 4> yyy = {0.0, 0.0, 3.0, 3.0};
    EXPECT_NEAR(jacobianDet(xxx, yyy, 0.1, 0.2), 1.5, 1e-12);
}

TEST(FeKernels, ApplyTensorAndDots)
{
    const Tensor ten = {2.0, 0.0, 0.0, 4.0};
    const Flux out = applyTensor(ten, Flux{3.0, 5.0});
    EXPECT_NEAR(out[0], 6.0, 1e-12);
    EXPECT_NEAR(out[1], 20.0, 1e-12);
    EXPECT_NEAR(dot2(Flux{1.0, 2.0}, Flux{3.0, 4.0}), 11.0, 1e-12);
    EXPECT_NEAR(dot4({1.0, 2.0, 3.0, 4.0}, {1.0, 1.0, 1.0, 1.0}), 10.0, 1e-12);
}
