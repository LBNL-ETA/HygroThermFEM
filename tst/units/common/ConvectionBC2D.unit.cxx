#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

TEST(TestConvectionBC2D, TestIntegrationPoints)
{
    SCOPED_TRACE("Begin Test: Convection boundary condition.");

    HygroThermFEM::MultiDomain multiDomain({.performThermal = true, .performMoisture = false});

    // Enter nodes. Arguments are: node number, x-coordinate, y-coordinate

    multiDomain.nodes().createNode({.index = 1, .x = 15, .y = 5});
    multiDomain.nodes().createNode({.index = 2, .x = 15, .y = 0});

    constexpr auto hc = 20.0;
    constexpr auto tAir = 255.15;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff{tAir, hc};

    auto aBc = HygroThermFEM::ConstantConvectionBC(multiDomain.nodes(), 1, 2, bcCoeff);

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
