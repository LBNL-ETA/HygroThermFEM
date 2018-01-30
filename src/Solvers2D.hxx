#pragma once

#include <vector>

#include "Elements2D.hxx"
#include "BoundaryConditions2D.hxx"

namespace MoisThermFEM {

  class ElementsLinear2D;
  class BoundaryConditions2D;

  /////////////////////////////////////////////////////////////////////////
  // ISolver2D
  /////////////////////////////////////////////////////////////////////////

  class ISolver2D {
  public:
    virtual ~ISolver2D() = default;
    ISolver2D( 
      ElementsLinear2D const & t_Elements,
      BoundaryConditions2D const & t_BCs );

    virtual void solve() = 0;

  protected:
    ElementsLinear2D m_Elements;
    BoundaryConditions2D m_BCs;
    bool m_Solved;
  };

  /////////////////////////////////////////////////////////////////////////
  // SteadyStateSolver2D
  /////////////////////////////////////////////////////////////////////////

  class SteadyStateSolver2D : public ISolver2D {
  public:
    SteadyStateSolver2D( 
      ElementsLinear2D const & t_Elements,
      BoundaryConditions2D const & t_BCs );

    std::vector< double > solution();

    void solve() override;

  private:
    std::vector< double > m_Solution;
  };

  /////////////////////////////////////////////////////////////////////////
  // TransientSolver2D
  /////////////////////////////////////////////////////////////////////////

  class TransientSolver2D : public ISolver2D {
  public:
    TransientSolver2D( 
      const ElementsLinear2D & t_Elements,
      const BoundaryConditions2D & t_BCs,
      const double DTime,
      const std::size_t NTimeSteps );

    std::vector< std::vector< double > > solution();

    void solve() override;

  private:
    std::vector< std::vector< double > > m_Solution;
    const double m_DTime;
    const std::size_t m_NSteps;
  };
  
}