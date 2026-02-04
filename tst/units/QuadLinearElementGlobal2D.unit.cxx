#include <memory>
#include <stdexcept>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

class TestQuadLinearElementGlobal2D : public testing::Test
{
protected:
    void SetUp() override
    {}

    void TearDown() override
    {}
};

TEST_F(TestQuadLinearElementGlobal2D, TestIntegrationPoint1)
{
    SCOPED_TRACE(
      "Begin Test: Quadrilateral linear element 2D in global coordinates - Gauss points.");

    HygroThermFEM::MultiDomain multiDomain(true, false);

    // Enter nodes. Arguments are: node number, x-coordinate, y-coordinate

    auto node1 = multiDomain.nodePool().createNode(1, 0, 0);
    auto node2 = multiDomain.nodePool().createNode(2, 5, 0);
    auto node3 = multiDomain.nodePool().createNode(3, 5, 5);
    auto node4 = multiDomain.nodePool().createNode(4, 0, 5);

    auto aElement = HygroThermFEM::QuadrilateralLinearGlobal2D(node1, node2, node3, node4);

    /////////////////////////////////////////////////////////
    ///    Integration Point 1
    /////////////////////////////////////////////////////////
    const size_t intPointIndex1 = 0;
    auto xg = aElement.xg(intPointIndex1);
    auto yg = aElement.yg(intPointIndex1);

    EXPECT_NEAR(1.056624327, xg, 1e-6);
    EXPECT_NEAR(1.056624327, yg, 1e-6);

    auto dPsiDx = aElement.DPsiDx(intPointIndex1);
    std::vector<double> correctDPsiDx{-0.157735027, 0.157735027, 0.042264973, -0.042264973};

    EXPECT_EQ(dPsiDx.size(), correctDPsiDx.size());

    for(auto i = 0u; i < correctDPsiDx.size(); ++i)
    {
        EXPECT_NEAR(correctDPsiDx[i], dPsiDx[i], 1e-6);
    }

    auto DPsiDy = aElement.DPsiDy(intPointIndex1);
    std::vector<double> correctDPsiDy{-0.157735027, -0.042264973, 0.042264973, 0.157735027};

    EXPECT_EQ(DPsiDy.size(), correctDPsiDy.size());

    for(auto i = 0u; i < correctDPsiDy.size(); ++i)
    {
        EXPECT_NEAR(correctDPsiDy[i], DPsiDy[i], 1e-6);
    }

    /////////////////////////////////////////////////////////
    ///    Integration Point 2
    /////////////////////////////////////////////////////////
    const size_t intPointIndex2 = 1;
    xg = aElement.xg(intPointIndex2);
    yg = aElement.yg(intPointIndex2);

    EXPECT_NEAR(3.943375673, xg, 1e-6);
    EXPECT_NEAR(1.056624327, yg, 1e-6);

    dPsiDx = aElement.DPsiDx(intPointIndex2);
    correctDPsiDx = {-0.157735027, 0.157735027, 0.042264973, -0.042264973};

    EXPECT_EQ(dPsiDx.size(), correctDPsiDx.size());

    for(auto i = 0u; i < correctDPsiDx.size(); ++i)
    {
        EXPECT_NEAR(correctDPsiDx[i], dPsiDx[i], 1e-6);
    }

    DPsiDy = aElement.DPsiDy(intPointIndex2);
    correctDPsiDy = {-0.042264973, -0.157735027, 0.157735027, 0.042264973};

    EXPECT_EQ(DPsiDy.size(), correctDPsiDy.size());

    for(auto i = 0u; i < correctDPsiDy.size(); ++i)
    {
        EXPECT_NEAR(correctDPsiDy[i], DPsiDy[i], 1e-6);
    }

    /////////////////////////////////////////////////////////
    ///    Integration Point 3
    /////////////////////////////////////////////////////////
    const size_t intPointIndex3 = 2;
    xg = aElement.xg(intPointIndex3);
    yg = aElement.yg(intPointIndex3);

    EXPECT_NEAR(3.943375673, xg, 1e-6);
    EXPECT_NEAR(3.943375673, yg, 1e-6);

    dPsiDx = aElement.DPsiDx(intPointIndex3);
    correctDPsiDx = {-0.042264973, 0.042264973, 0.157735027, -0.157735027};

    EXPECT_EQ(dPsiDx.size(), correctDPsiDx.size());

    for(auto i = 0u; i < correctDPsiDx.size(); ++i)
    {
        EXPECT_NEAR(correctDPsiDx[i], dPsiDx[i], 1e-6);
    }

    DPsiDy = aElement.DPsiDy(intPointIndex3);
    correctDPsiDy = {-0.042264973, -0.157735027, 0.157735027, 0.042264973};

    EXPECT_EQ(DPsiDy.size(), correctDPsiDy.size());

    for(auto i = 0u; i < correctDPsiDy.size(); ++i)
    {
        EXPECT_NEAR(correctDPsiDy[i], DPsiDy[i], 1e-6);
    }

    /////////////////////////////////////////////////////////
    ///    Integration Point 4
    /////////////////////////////////////////////////////////
    const size_t intPointIndex4 = 3;
    xg = aElement.xg(intPointIndex4);
    yg = aElement.yg(intPointIndex4);

    EXPECT_NEAR(1.056624327, xg, 1e-6);
    EXPECT_NEAR(3.943375673, yg, 1e-6);

    dPsiDx = aElement.DPsiDx(intPointIndex4);
    correctDPsiDx = {-0.042264973, 0.042264973, 0.157735027, -0.157735027};

    EXPECT_EQ(dPsiDx.size(), correctDPsiDx.size());

    for(auto i = 0u; i < correctDPsiDx.size(); ++i)
    {
        EXPECT_NEAR(correctDPsiDx[i], dPsiDx[i], 1e-6);
    }

    DPsiDy = aElement.DPsiDy(intPointIndex4);
    correctDPsiDy = {-0.157735027, -0.042264973, 0.042264973, 0.157735027};

    EXPECT_EQ(DPsiDy.size(), correctDPsiDy.size());

    for(auto i = 0u; i < correctDPsiDy.size(); ++i)
    {
        EXPECT_NEAR(correctDPsiDy[i], DPsiDy[i], 1e-6);
    }
}
