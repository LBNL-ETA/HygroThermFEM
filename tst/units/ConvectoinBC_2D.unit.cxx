#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

class TestBoundaryConditions2D_test1 : public testing::Test
{
protected:
    void SetUp() override
    {}

    void TearDown() override
    {}
};

TEST_F(TestBoundaryConditions2D_test1, TestIntegrationPoints)
{
    SCOPED_TRACE("Begin Test: Convection boundary condition integral test.");

    HygroThermFEM::MultiDomain multiDomain(true, false);

    // Enter nodes. Arguments are: node number, x-coordinate, y-coordinate
    multiDomain.nodePool().createNode(1, 15, 5);
    multiDomain.nodePool().createNode(2, 15, 0);
    multiDomain.nodePool().createNode(3, 5, 5);
    multiDomain.nodePool().createNode(4, 5, 0);
    multiDomain.nodePool().createNode(5, 0, 5);
    multiDomain.nodePool().createNode(6, 0, 0);

    auto const hc1 = 20.0;
    auto const Tair1 = 255.15;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff1{Tair1, hc1};

    auto const hc2 = 2.4;
    auto const Tair2 = 294.15;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff2{Tair2, hc2};

    HygroThermFEM::BoundaryConditions2D BCs;
    BCs.assignBC(fem::make_unique<HygroThermFEM::ConstantConvectionBC>(multiDomain.nodePool(), 1, 2, bcCoeff1));
    BCs.assignBC(fem::make_unique<HygroThermFEM::ConstantConvectionBC>(multiDomain.nodePool(), 6, 5, bcCoeff2));

    auto maxNodeIndex = multiDomain.nodePool().maxIndex();
    auto HMat = BCs.HMatrix(maxNodeIndex);

    std::vector<std::vector<double>> correctH{{33.33333333, 16.66666667, 0, 0, 0, 0},
                                              {16.66666667, 33.33333333, 0, 0, 0, 0},
                                              {0, 0, 0, 0, 0, 0},
                                              {0, 0, 0, 0, 0, 0},
                                              {0, 0, 0, 0, 4, 2},
                                              {0, 0, 0, 0, 2, 4}};
    EXPECT_EQ(correctH.size(), HMat.size());

    for(auto idx = 0u; idx < correctH.size(); ++idx)
    {
        for(auto jdx = 0u; jdx < correctH.size(); ++jdx)
        {
            EXPECT_NEAR(correctH[idx][jdx], HMat(idx, jdx), 1e-6);
        }
    }

    auto vecR = BCs.RVector(maxNodeIndex);

    std::vector<double> correctR{12757.5, 12757.5, 0, 0, 1764.9, 1764.9};

    EXPECT_EQ(vecR.size(), correctR.size());

    for(auto idx = 0u; idx < correctR.size(); ++idx)
    {
        EXPECT_NEAR(correctR[idx], vecR[idx], 1e-6);
    }
}
