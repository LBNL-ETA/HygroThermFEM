#ifndef LINEARSOLVER_H
#define LINEARSOLVER_H

#include <memory>
#include <vector>

#include "SquareMatrix.hxx"

namespace FenestrationCommon {

  class CLinearSolver {
  public:
    CLinearSolver();

    std::vector< double > solveSystem( 
      FenestrationCommon::SquareMatrix< double > t_MatrixA, 
      std::vector< double > & t_VectorB ) const;

  private:
    std::vector< double > checkSingularity( SquareMatrix< double > & t_MatrixA ) const;
    std::vector< size_t > makeUpperTriangular( SquareMatrix< double > & t_MatrixA ) const;

  };
}

#endif