#pragma once

#include <vector>
#include <memory>

#include "LocalShapeFunctions.hxx"

namespace MoisThermFEM {
  
  ////////////////////////////////////////////////////////////////////////////
  //   QuadrilateralLinearLocal2D
  ////////////////////////////////////////////////////////////////////////////

  // Represents quadrilateral linear element in local coordinate system (-1, 1)
  class QuadrilateralLinearLocal2D {
  public:
    static QuadrilateralLinearLocal2D& Instance();

    double Psi( size_t const IntegrationPointIndex, size_t const Index );
    double PsiDKsi( size_t const IntegrationPointIndex, size_t const Index );
    double PsiDEta( size_t const IntegrationPointIndex, size_t const Index );

    std::vector< double > VPsi( size_t const IntegrationPointIndex ) const;
    std::vector< double > VPsiDKsi( size_t const IntegrationPointIndex ) const;
    std::vector< double > VPsiDEta( size_t const IntegrationPointIndex ) const;

  private:
    QuadrilateralLinearLocal2D();
    ~QuadrilateralLinearLocal2D();

    // Shape functions for every integration point (at local coordinate system)
    std::vector< std::unique_ptr< ILocalShapeFunctions2DQuadrilateral > > m_Ksi;

    // Class defines shape functions and its derivatives in local coordinate system
    // for quadrilateral linear element.
    class QuadLinearLocalShapeFunctions2D : public ILocalShapeFunctions2DQuadrilateral {
    public:
      QuadLinearLocalShapeFunctions2D( LocalPoint2D const t_Point );

    };

  };
  
}