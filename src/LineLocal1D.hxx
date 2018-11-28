#pragma once

#include <vector>
#include <memory>

#include "LocalShapeFunctions.hxx"

namespace MoisThermFEM {
  
  class LineLinearLocal1D {
  public:
    static LineLinearLocal1D& Instance();

    double Psi( size_t IntegrationPointIndex, size_t Index );

  private:
    LineLinearLocal1D();

    std::vector< std::unique_ptr< ILocalShapeFunctions1DLine > > m_Ksi;
		// std::vector< std::reference_wrapper< ILocalShapeFunctions1DLine > > m_Ksi;

    class LineLinearLocalShapeFunctions1D : public ILocalShapeFunctions1DLine {
    public:
      explicit LineLinearLocalShapeFunctions1D( LocalPoint1D const & t_Point );
    };
  };

}