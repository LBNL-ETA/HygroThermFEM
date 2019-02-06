#include <gtest/gtest.h>

#include "MoisThermFEM2D.hxx"

using MoisThermFEM::NodePool;


class TestConvectionBC2D : public testing::Test
{
protected:
    void SetUp() override
    {}

    void TearDown() override
    {
        NodePool::Instance().clear();
    }
};

TEST_F(TestConvectionBC2D, TestIntegrationPoints)
{
    SCOPED_TRACE("Begin Test: Convection boundary condition.");

    // Enter nodes. Arguments are: node number, x-coordinate, y-coordinate

    NodePool::Instance().createNode(1, 15, 5);
    NodePool::Instance().createNode(2, 15, 0);

    const auto hc = 20.0;
    const auto tAir = 255.15;

    auto aBc = MoisThermFEM::ConstantConvectionBC(1, 2, tAir, hc);

    auto h = aBc.H_Matrix();

    std::vector<std::vector<double>> correctH{{33.33333333, 16.66666667},
                                              {16.66666667, 33.33333333}};
    EXPECT_EQ(correctH.size(), h.size());

    for(auto i = 0u; i < correctH.size(); ++i)
    {
        for(auto j = 0u; j < correctH.size(); ++j)
        {
            EXPECT_NEAR(correctH[i][j], h(i, j), 1e-6);
        }
    }

    auto R = aBc.R_Vector();

    std::vector<double> correctR{12757.5, 12757.5};

    EXPECT_EQ(R.size(), correctR.size());

    for(auto i = 0u; i < correctR.size(); ++i)
    {
        EXPECT_NEAR(correctR[i], R[i], 1e-6);
    }
}
