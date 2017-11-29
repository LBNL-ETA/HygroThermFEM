#include <stdexcept>
#include <assert.h>

#include "LinearSolver.hxx"
#include "SquareMatrix.hxx"

#include <math.h>

namespace FenestrationCommon {

  CLinearSolver::CLinearSolver() {
  
  }

  std::vector< double > CLinearSolver::checkSingularity( CSquareMatrix& t_MatrixA ) const {
    auto size = t_MatrixA.getSize();
    std::vector< double > vv;

    for( size_t i = 0 ; i < size; ++i ) {
      double aamax = 0;
      for( size_t j = 0 ; j < size; ++j ) {
        auto absCellValue = fabs( t_MatrixA[ i ][ j ] );
        if ( absCellValue > aamax ) {
          aamax = absCellValue;
        }
      }
      if( aamax == 0 ) {
        assert( aamax != 0 );
      }
      vv.push_back( 1 / aamax );
    }

    return vv;
  }

  std::vector< size_t > CLinearSolver::makeUpperTriangular( CSquareMatrix& t_MatrixA ) const {
    auto TINY = 1e-20;

    auto size = int( t_MatrixA.getSize() );
    std::vector< size_t > index( size );

    auto vv = checkSingularity( t_MatrixA );

    auto d = 1;

    for( auto j = 0; j < size; ++j ) {

      for( auto i = 0; i <= j - 1; ++i ) {
        auto sum = t_MatrixA[ i ][ j ];
        for( auto k = 0; k <= i - 1; ++k ) {
          sum = sum - t_MatrixA[ i ][ k ] * t_MatrixA[ k ][ j ];
        }
        t_MatrixA[ i ][ j ] = sum;
      }

      double aamax = 0;
      auto imax = 0;

      for( auto i = j; i < size; ++i ) {
        auto sum = t_MatrixA[ i ][ j ];
        for( auto k = 0; k <= j - 1; ++k ) {
          sum = sum - t_MatrixA[ i ][ k ] * t_MatrixA[ k ][ j ];
        }
        t_MatrixA[ i ][ j ] = sum;
        auto dum = vv[i] * fabs( sum );
        if ( dum >= aamax ) {
          imax = i;
          aamax = dum;
        }
      }

      if ( j != imax ) {
        for ( auto k = 0; k < size; ++k ) {
          auto dum = t_MatrixA[ imax ][ k ];
          t_MatrixA[ imax ][ k ] = t_MatrixA[ j ][ k ];
          t_MatrixA[ j ][ k ] = dum;
        } // k
        d = -d;
        vv[ imax ] = vv[ j ];
      }
      index[ j ] = imax;
      if ( t_MatrixA[ j ][ j ] == 0.0 ) {
        t_MatrixA[ j ][ j ] = TINY;
      }
      if ( j != ( size - 1 ) ) {
        auto dum = 1.0 / t_MatrixA[ j ][ j ];
        for( auto i = j + 1; i < size; ++i ) {
          t_MatrixA[ i ][ j ] = t_MatrixA[ i ][ j ] * dum;
        } // i
      }

    }

    return index;
  }

  std::vector< double > CLinearSolver::solveSystem(
    CSquareMatrix t_MatrixA, 
    std::vector< double >& t_VectorB ) const {

    if ( t_MatrixA.getSize() != t_VectorB.size() ) {
      std::runtime_error( "Matrix and vector for system of linear equations are not same size." );
    }

    auto index = makeUpperTriangular( t_MatrixA );

    auto size = int( t_MatrixA.getSize() );

    auto ii = -1;
    for( auto i = 0; i < size; ++i ) {
      auto ll = index[ i ];
      auto sum = t_VectorB[ ll ];
      t_VectorB[ ll ] = t_VectorB[ i ];
      if( ii != -1 ) {
        for( auto j = ii; j <= i - 1; ++j ) {
          sum -= t_MatrixA[ i ][ j ] * t_VectorB[ j ];
        } // j
      } else if ( sum != 0.0 ) {
        ii = int( i );
      }
      t_VectorB[ i ] = sum;
    } // i

    for( auto i = ( size - 1 ); i >= 0; --i ) {
      auto sum = t_VectorB[ i ];
      for( auto j = i + 1; j < size; ++j ) {
        sum -= t_MatrixA[ i ][ j ] * t_VectorB[ j ];
      } // j
      t_VectorB[ i ] = sum / t_MatrixA[ i ][ i ];
    } // i

    return t_VectorB;

  }

}