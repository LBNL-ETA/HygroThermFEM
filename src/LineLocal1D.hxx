#pragma once

#include <vector>
#include <memory>

#include "LocalShapeFunctions.hxx"

namespace Conrad {
  
  class LineLinearLocal1D {
  public:
    static LineLinearLocal1D& Instance();

    double Psi( size_t const IntegrationPointIndex, size_t const Index );

  private:
    LineLinearLocal1D();
    ~LineLinearLocal1D();

    std::vector< std::shared_ptr< ILocalShapeFunctions1DLine > > m_Ksi;

    class LineLinearLocalShapeFunctions1D : public ILocalShapeFunctions1DLine {
    public:
      explicit LineLinearLocalShapeFunctions1D( LocalPoint1D const & t_Point );
    };
  };

}