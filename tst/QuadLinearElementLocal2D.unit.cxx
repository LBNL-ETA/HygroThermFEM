#include <memory>
#include <stdexcept>
#include <gtest/gtest.h>

#include <memory>

#include "MoisThermFEM2D.hxx"

using namespace MoisThermFEM;

class TestQuadLinearElementLocal2D : public testing::Test {

protected:
    void
    SetUp() override {
    }

};

TEST_F( TestQuadLinearElementLocal2D, TestIntegrationPoint1 ) {
    SCOPED_TRACE( "Begin Test: Quadrilateral linear element 2D in local coordinates." );

    auto & aElement = QuadrilateralLinearLocal2D::Instance();

    // Integration point 1
    const auto IntegrationPoint = 0;

    // Psi
    std::vector< double > correctPsi = { 0.622008468, 0.166666667, 0.044658199, 0.166666667 };
    auto psi = aElement.VPsi( IntegrationPoint );

    EXPECT_EQ( correctPsi.size(), psi.size() );

    for ( size_t i = 0u; i < correctPsi.size(); ++i ) {
        EXPECT_NEAR( correctPsi[ i ], psi[ i ], 1e-6 );
    }

    // PsiDKsi
    std::vector< double > correctPsiDKsi = { -0.394337567, 0.394337567, 0.105662433, -0.105662433 };
    auto psiDKsi = aElement.VPsiDKsi( IntegrationPoint );

    EXPECT_EQ( correctPsiDKsi.size(), psiDKsi.size() );

    for ( size_t i = 0u; i < correctPsiDKsi.size(); ++i ) {
        EXPECT_NEAR( correctPsiDKsi[ i ], psiDKsi[ i ], 1e-6 );
    }

    // PsiDEta
    std::vector< double > correctPsiDEta = { -0.394337567, -0.105662433, 0.105662433, 0.394337567 };
    auto psiDEta = aElement.VPsiDEta( IntegrationPoint );

    EXPECT_EQ( correctPsiDEta.size(), psiDEta.size() );

    for ( size_t i = 0u; i < correctPsiDEta.size(); ++i ) {
        EXPECT_NEAR( correctPsiDEta[ i ], psiDEta[ i ], 1e-6 );
    }
}

TEST_F( TestQuadLinearElementLocal2D, TestIntegrationPoint2 ) {
    SCOPED_TRACE( "Begin Test: Quadrilateral linear element 2D in local coordinates." );

    auto & aElement = QuadrilateralLinearLocal2D::Instance();

    // Integration point 2
    const auto IntegrationPoint = 1;

    // Psi
    std::vector< double > correctPsi = { 0.166666667, 0.622008468, 0.166666667, 0.044658199 };
    auto psi = aElement.VPsi( IntegrationPoint );

    EXPECT_EQ( correctPsi.size(), psi.size() );

    for ( auto i = 0u; i < correctPsi.size(); ++i ) {
        EXPECT_NEAR( correctPsi[ i ], psi[ i ], 1e-6 );
    }

    // PsiDKsi
    std::vector< double > correctPsiDKsi = { -0.394337567, 0.394337567, 0.105662433, -0.105662433 };
    auto psiDKsi = aElement.VPsiDKsi( IntegrationPoint );

    EXPECT_EQ( correctPsiDKsi.size(), psiDKsi.size() );

    for ( auto i = 0u; i < correctPsiDKsi.size(); ++i ) {
        EXPECT_NEAR( correctPsiDKsi[ i ], psiDKsi[ i ], 1e-6 );
    }

    // PsiDEta
    std::vector< double > correctPsiDEta = { -0.105662433, -0.394337567, 0.394337567, 0.105662433 };
    auto PsiDEta = aElement.VPsiDEta( IntegrationPoint );

    EXPECT_EQ( correctPsiDEta.size(), PsiDEta.size() );

    for ( auto i = 0u; i < correctPsiDEta.size(); ++i ) {
        EXPECT_NEAR( correctPsiDEta[ i ], PsiDEta[ i ], 1e-6 );
    }
}

TEST_F( TestQuadLinearElementLocal2D, TestIntegrationPoint3 ) {
    SCOPED_TRACE( "Begin Test: Quadrilateral linear element 2D in local coordinates." );

    auto & aElement = QuadrilateralLinearLocal2D::Instance();

    // Integration point 3
    const auto IntegrationPoint = 2;

    // Psi
    std::vector< double > correctPsi = { 0.044658199, 0.166666667, 0.622008468, 0.166666667 };
    auto psi = aElement.VPsi( IntegrationPoint );

    EXPECT_EQ( correctPsi.size(), psi.size() );

    for ( auto i = 0u; i < correctPsi.size(); ++i ) {
        EXPECT_NEAR( correctPsi[ i ], psi[ i ], 1e-6 );
    }

    // PsiDKsi
    std::vector< double > correctPsiDKsi = { -0.105662433, 0.105662433, 0.394337567, -0.394337567 };
    auto psiDKsi = aElement.VPsiDKsi( IntegrationPoint );

    EXPECT_EQ( correctPsiDKsi.size(), psiDKsi.size() );

    for ( auto i = 0u; i < correctPsiDKsi.size(); ++i ) {
        EXPECT_NEAR( correctPsiDKsi[ i ], psiDKsi[ i ], 1e-6 );
    }

    // PsiDEta
    std::vector< double > correctPsiDEta = { -0.105662433, -0.394337567, 0.394337567, 0.105662433 };
    auto psiDEta = aElement.VPsiDEta( IntegrationPoint );

    EXPECT_EQ( correctPsiDEta.size(), psiDEta.size() );

    for ( auto i = 0u; i < correctPsiDEta.size(); ++i ) {
        EXPECT_NEAR( correctPsiDEta[ i ], psiDEta[ i ], 1e-6 );
    }
}

TEST_F( TestQuadLinearElementLocal2D, TestIntegrationPoint4 ) {
    SCOPED_TRACE( "Begin Test: Quadrilateral linear element 2D in local coordinates." );

    auto & aElement = QuadrilateralLinearLocal2D::Instance();

    // Integration point 4
    const auto integrationPoint = 3;

    // Psi
    std::vector< double > correctPsi = { 0.166666667, 0.044658199, 0.166666667, 0.622008468 };
    auto psi = aElement.VPsi( integrationPoint );

    EXPECT_EQ( correctPsi.size(), psi.size() );

    for ( auto i = 0u; i < correctPsi.size(); ++i ) {
        EXPECT_NEAR( correctPsi[ i ], psi[ i ], 1e-6 );
    }

    // PsiDKsi
    std::vector< double > correctPsiDKsi = { -0.105662433, 0.105662433, 0.394337567, -0.394337567 };
    auto psiDKsi = aElement.VPsiDKsi( integrationPoint );

    EXPECT_EQ( correctPsiDKsi.size(), psiDKsi.size() );

    for ( auto i = 0u; i < correctPsiDKsi.size(); ++i ) {
        EXPECT_NEAR( correctPsiDKsi[ i ], psiDKsi[ i ], 1e-6 );
    }

    // PsiDEta
    std::vector< double > correctPsiDEta = { -0.394337567, -0.105662433, 0.105662433, 0.394337567 };
    auto psiDEta = aElement.VPsiDEta( integrationPoint );

    EXPECT_EQ( correctPsiDEta.size(), psiDEta.size() );

    for ( auto i = 0u; i < correctPsiDEta.size(); ++i ) {
        EXPECT_NEAR( correctPsiDEta[ i ], psiDEta[ i ], 1e-6 );
    }
}
