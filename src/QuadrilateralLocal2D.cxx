#include "QuadrilateralLocal2D.hxx"
#include "Node2D.hxx"
#include "IntegrationPoints.hxx"

namespace MoisThermFEM {

  ////////////////////////////////////////////////////////////////////////////
  //   QuadElement2DLinearLocal::QuadLinearPoint2D
  ////////////////////////////////////////////////////////////////////////////

  QuadrilateralLinearLocal2D::QuadLinearLocalShapeFunctions2D::QuadLinearLocalShapeFunctions2D(
    LocalPoint2D const t_Point ) :
    ILocalShapeFunctions2DQuadrilateral( t_Point ) {

    // Form shape functions
    m_Psi.push_back( 0.25 * ( 1 - t_Point.ksi ) * ( 1 - t_Point.eta ) );
    m_Psi.push_back( 0.25 * ( 1 + t_Point.ksi ) * ( 1 - t_Point.eta ) );
    m_Psi.push_back( 0.25 * ( 1 + t_Point.ksi ) * ( 1 + t_Point.eta ) );
    m_Psi.push_back( 0.25 * ( 1 - t_Point.ksi ) * ( 1 + t_Point.eta ) );

    // Form derivatives per ksi
    m_PsiDKsi.push_back( -0.25 * ( 1 - t_Point.eta ) );
    m_PsiDKsi.push_back( 0.25 * ( 1 - t_Point.eta ) );
    m_PsiDKsi.push_back( 0.25 * ( 1 + t_Point.eta ) );
    m_PsiDKsi.push_back( -0.25 * ( 1 + t_Point.eta ) );

    // Form derivatives per eta
    m_PsiDEta.push_back( -0.25 * ( 1 - t_Point.ksi ) );
    m_PsiDEta.push_back( -0.25 * ( 1 + t_Point.ksi ) );
    m_PsiDEta.push_back( 0.25 * ( 1 + t_Point.ksi ) );
    m_PsiDEta.push_back( 0.25 * ( 1 - t_Point.ksi ) );
  }

  ////////////////////////////////////////////////////////////////////////////
  //   QuadrilateralLinearLocal2D
  ////////////////////////////////////////////////////////////////////////////

  QuadrilateralLinearLocal2D& QuadrilateralLinearLocal2D::Instance() {
    //Thread safe in C++11
    static QuadrilateralLinearLocal2D m_Instance;
    return m_Instance;
  }

  QuadrilateralLinearLocal2D::QuadrilateralLinearLocal2D() {
    auto aPoints = IntegrationPoints2D::Instance().getPoints2D();
    for( auto point : aPoints ) {
      m_Ksi.push_back( std::unique_ptr< QuadLinearLocalShapeFunctions2D >( new QuadLinearLocalShapeFunctions2D( point ) ) );
    }
  }

  QuadrilateralLinearLocal2D::~QuadrilateralLinearLocal2D() {

  }

  double QuadrilateralLinearLocal2D::Psi( size_t const IntegrationPointIndex, size_t const Index ) {
    if( IntegrationPointIndex >= m_Ksi.size() ) {
      throw std::runtime_error( "Integration point index out of range. Rouinte QuadElement2DLinearLocal::Psi." );
    }

    return m_Ksi[ IntegrationPointIndex ]->Psi( Index );
  }

  double QuadrilateralLinearLocal2D::PsiDKsi( size_t const IntegrationPointIndex, size_t const Index ) {
    if( IntegrationPointIndex >= m_Ksi.size() ) {
      throw std::runtime_error( "Integration point index out of range. Rouinte QuadElement2DLinearLocal::PsiDKsi." );
    }

    return m_Ksi[ IntegrationPointIndex ]->PsiDKsi( Index );
  }

  double QuadrilateralLinearLocal2D::PsiDEta( size_t const IntegrationPointIndex, size_t const Index ) {
    if( IntegrationPointIndex >= m_Ksi.size() ) {
      throw std::runtime_error( "Integration point index out of range. Rouinte QuadElement2DLinearLocal::PsiDEta." );
    }

    return m_Ksi[ IntegrationPointIndex ]->PsiDEta( Index );
  }

  std::vector< double > QuadrilateralLinearLocal2D::VPsi( size_t const IntegrationPointIndex ) const {
    if( IntegrationPointIndex >= m_Ksi.size() ) {
      throw std::runtime_error( "Integration point index out of range. Rouinte QuadElement2DLinearLocal::VPsi." );
    }

    return m_Ksi[ IntegrationPointIndex ]->VPsi();
  }

  std::vector< double > QuadrilateralLinearLocal2D::VPsiDKsi( size_t const IntegrationPointIndex ) const {
    if( IntegrationPointIndex >= m_Ksi.size() ) {
      throw std::runtime_error( "Integration point index out of range. Rouinte QuadElement2DLinearLocal::VPsiDKsi." );
    }

    return m_Ksi[ IntegrationPointIndex ]->VPsiDKsi();
  }

  std::vector< double > QuadrilateralLinearLocal2D::VPsiDEta( size_t const IntegrationPointIndex ) const {
    if( IntegrationPointIndex >= m_Ksi.size() ) {
      throw std::runtime_error( "Integration point index out of range. Rouinte QuadElement2DLinearLocal::VPsiDEta." );
    }

    return m_Ksi[ IntegrationPointIndex ]->VPsiDEta();
  }

}