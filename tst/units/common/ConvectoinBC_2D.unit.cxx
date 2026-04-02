#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

TEST(TestBoundaryConditions2D_test1, TestIntegrationPoints)
{
    SCOPED_TRACE("Begin Test: Convection boundary condition integral test.");

    HygroThermFEM::MultiDomain multiDomain({.performThermal = true, .performMoisture = false});

    // Enter nodes. Arguments are: node number, x-coordinate, y-coordinate
    multiDomain.nodes().createNode({.index = 1, .x = 15, .y = 5});
    multiDomain.nodes().createNode({.index = 2, .x = 15, .y = 0});
    multiDomain.nodes().createNode({.index = 3, .x = 5, .y = 5});
    multiDomain.nodes().createNode({.index = 4, .x = 5, .y = 0});
    multiDomain.nodes().createNode({.index = 5, .x = 0, .y = 5});
    multiDomain.nodes().createNode({.index = 6, .x = 0, .y = 0});

    auto const hc1 = 20.0;
    auto const Tair1 = 255.15;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff1{Tair1, hc1};

    auto const hc2 = 2.4;
    auto const Tair2 = 294.15;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff2{Tair2, hc2};

    HygroThermFEM::BoundaryConditions2D BCs;
    BCs.assignBC(fem::make_unique<HygroThermFEM::ConstantConvectionBC>(multiDomain.nodes(), 1, 2, bcCoeff1));
    BCs.assignBC(fem::make_unique<HygroThermFEM::ConstantConvectionBC>(multiDomain.nodes(), 6, 5, bcCoeff2));

    auto maxNodeIndex = multiDomain.nodes().maxIndex();
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
