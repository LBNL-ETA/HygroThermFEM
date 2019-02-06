#include <gtest/gtest.h>

#include "MoisThermFEM2D.hxx"

using HygroThermFEM::NodePool;
using HygroThermFEM::MaterialPool;

class TestBoundaryConditions2D_test1 : public testing::Test
{
protected:
    void SetUp() override
    {}

    void TearDown() override
    {
        NodePool::Instance().clear();
        MaterialPool::Instance().clear();
    }
};

TEST_F(TestBoundaryConditions2D_test1, TestIntegrationPoints)
{
    SCOPED_TRACE("Begin Test: Convection boundary condition integral test.");

    // Enter nodes. Arguments are: node number, x-coordinate, y-coordinate
    NodePool::Instance().createNode(1, 15, 5);
    NodePool::Instance().createNode(2, 15, 0);
    NodePool::Instance().createNode(3, 5, 5);
    NodePool::Instance().createNode(4, 5, 0);
    NodePool::Instance().createNode(5, 0, 5);
    NodePool::Instance().createNode(6, 0, 0);

    auto const hc1 = 20.0;
    auto const Tair1 = 255.15;

    auto const hc2 = 2.4;
    auto const Tair2 = 294.15;

    HygroThermFEM::BoundaryConditions2D BCs;
    BCs.assignBC(fem::make_unique<HygroThermFEM::ConstantConvectionBC>(1, 2, Tair1, hc1));
    BCs.assignBC(fem::make_unique<HygroThermFEM::ConstantConvectionBC>(6, 5, Tair2, hc2));

    auto H = BCs.HMatrix();

    std::vector<std::vector<double>> correctH{{33.33333333, 16.66666667, 0, 0, 0, 0},
                                              {16.66666667, 33.33333333, 0, 0, 0, 0},
                                              {0, 0, 0, 0, 0, 0},
                                              {0, 0, 0, 0, 0, 0},
                                              {0, 0, 0, 0, 4, 2},
                                              {0, 0, 0, 0, 2, 4}};
    EXPECT_EQ(correctH.size(), H.size());

    for(auto i = 0u; i < correctH.size(); ++i)
    {
        for(auto j = 0u; j < correctH.size(); ++j)
        {
            EXPECT_NEAR(correctH[i][j], H(i, j), 1e-6);
        }
    }

    auto R = BCs.RVector();

    std::vector<double> correctR{12757.5, 12757.5, 0, 0, 1764.9, 1764.9};

    EXPECT_EQ(R.size(), correctR.size());

    for(auto i = 0u; i < correctR.size(); ++i)
    {
        EXPECT_NEAR(correctR[i], R[i], 1e-6);
    }
}
