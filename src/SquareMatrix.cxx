#include <assert.h>
#include <stdexcept>
#include <thread>

#include "SquareMatrix.hxx"

namespace FenestrationCommon {
  CSquareMatrix::CSquareMatrix( size_t const aSize ) {
    m_Size = aSize;
    m_Matrix.resize( m_Size );
    for( unsigned i = 0; i < m_Size; ++i )
      m_Matrix[ i ].resize( m_Size );
  }

  size_t CSquareMatrix::getSize() const {
    return m_Size;
  }

  void CSquareMatrix::setZeros() {
    for( unsigned i = 0; i < m_Size; ++i ) {
      for( unsigned j = 0; j < m_Size; ++j ) {
        m_Matrix[ i ][ j ] = 0;
      }
    }
  }

  void CSquareMatrix::setIdentity() {
    setZeros();
    for( unsigned i = 0; i < m_Size; ++i ) {
      m_Matrix[ i ][ i ] = 1;
    }
  }

  void CSquareMatrix::setDiagonal( std::vector< double > const & t_Values ) {
    if( t_Values.size() != m_Size ) {
      throw std::runtime_error("Supplied vector size mismatch matrix size");
    }
    setZeros();
    for( unsigned i = 0; i < m_Size; ++i ) {
      m_Matrix[ i ][ i ] = t_Values[ i ];
    }
  }

  std::vector< double >& CSquareMatrix::operator[]( size_t const index ) {
    assert( index < m_Size );
    return m_Matrix[ index ];
  }

  CSquareMatrix CSquareMatrix::add( CSquareMatrix const & t_Matrix ) const {
    if( m_Size != t_Matrix.m_Size ) {
      throw std::runtime_error("Matrices must be identical in size.");
    }

		CSquareMatrix aMatrix{ m_Size };
    for( unsigned i = 0; i < m_Size; ++i ) {
      for( unsigned j = 0; j < t_Matrix.m_Size; ++j ) {
        aMatrix[ i ][ j ] = m_Matrix[ i ][ j ] + t_Matrix.m_Matrix[ i ][ j ];
      }
    }

    return aMatrix;
  }

  CSquareMatrix CSquareMatrix::addDiagonal( std::vector< double > const & t_Vector ) {
    if( m_Size != t_Vector.size() ) {
      std::runtime_error( "Matrix and vector have different sizes." );
    }

		CSquareMatrix aMatrix{ m_Size };
    for( unsigned i = 0; i < m_Size; ++i ) {
      for( unsigned j = 0; j < m_Size; ++j ) {
        aMatrix[ i ][ j ] = m_Matrix[ i ][ j ];
      }
      aMatrix[ i ][ i ] += t_Vector[ i ];
    }
    return aMatrix;
  }

  CSquareMatrix CSquareMatrix::sub( const CSquareMatrix& t_Matrix ) const {
    if( m_Size != t_Matrix.m_Size ) {
      throw std::runtime_error("Matrices must be identical in size.");
    }

		CSquareMatrix aMatrix{ m_Size };
    for( size_t i = 0; i < m_Size; ++i ) {
      for( size_t j = 0; j < t_Matrix.m_Size; ++j ) {
        aMatrix[ i ][ j ] = m_Matrix[ i ][ j ] - t_Matrix.m_Matrix[ i ][ j ];
      }
    }

    return aMatrix;
  }

  CSquareMatrix CSquareMatrix::mult( const CSquareMatrix& t_Matrix ) const {
    if( m_Size != t_Matrix.m_Size ) {
      throw std::runtime_error("Matrices must be identical in size.");
    }

		CSquareMatrix aMatrix{ m_Size };

    for( unsigned i = 0; i < m_Size; ++i ) {
      for( unsigned k = 0; k < m_Size; ++k ) {
        for( unsigned j = 0; j < t_Matrix.m_Size; ++j ) {
          aMatrix[ i ][ j ] = aMatrix[ i ][ j ] + m_Matrix[ i ][ k ] * t_Matrix.m_Matrix[ k ][ j ];
        }
      }
    }

    return aMatrix;

  }

  std::vector< double > CSquareMatrix::multMxV( std::vector< double > const & t_Vector ) const {
    if( m_Size != t_Vector.size() ) {
      throw std::runtime_error("Matrix and vector does not have same number of rows and columns."
        " It is not possible to perform multiplication.");
    }

		std::vector< double > aResult( m_Size );

    for( unsigned i = 0; i < m_Size; ++i ) {
      for( unsigned j = 0; j < m_Size; ++j ) {
        aResult[ i ] = aResult[ i ] + m_Matrix[ i ][ j ] * t_Vector[ j ];
      }
    }

    return aResult;
  }

  std::vector< double > CSquareMatrix::multVxM( const std::vector< double >& t_Vector ) const {
    if( m_Size != t_Vector.size() ) {
      throw std::runtime_error("Matrix and vector do not have same number of rows and columns."
                          " It is not possible to perform multiplication.");
    }

    std::vector< double > aResult( m_Size );

    for( unsigned i = 0; i < m_Size; ++i ) {
      for( unsigned j = 0; j < m_Size; ++j ) {
        aResult[ i ] = aResult[ i ] + m_Matrix[ j ][ i ] * t_Vector[ j ];
      }
    }

    return aResult;
  }

  void CSquareMatrix::copyFrom( const CSquareMatrix& t_Matrix ) {
    if( m_Size != t_Matrix.m_Size ) {
      throw std::runtime_error("Matrices must be identical in size");
    }
    for( unsigned i = 0; i < m_Size; ++i ) {
      m_Matrix[ i ] = t_Matrix.m_Matrix[ i ];
    }
  }

  CSquareMatrix CSquareMatrix::inverse() {
    // return LU decomposed matrix of current matrix
    auto aLU = LU();

    // find the inverse
		CSquareMatrix inverse{ m_Size };
    std::vector< double > d( m_Size );
    std::vector< double > y( m_Size );

	  const auto size = int( m_Size - 1 );

    for( auto m = 0; m <= size; ++m ) {
      std::fill( d.begin(), d.end(), 0 );
      std::fill( y.begin(), y.end(), 0 );
      d[ m ] = 1;
      for( auto i = 0; i <= size; ++i ) { 
        double x = 0;
        for( auto j = 0; j <= i - 1; ++j ) {
          x = x + aLU[ i ][ j ] * y[ j ];
        }
         y[ i ] = ( d[ i ] - x );
      }
  
      for( auto i = size; i >= 0; --i ) { 
        auto x = 0.0;
        for( auto j = i + 1; j <= size; ++j ) { 
          x = x + aLU[ i ][ j ] * inverse[ j ][ m ];
        }
         inverse[ i ][ m ] = ( y[ i ] - x ) / aLU[ i ][ i ];
      }
    }

    return inverse;
  }

    CSquareMatrix CSquareMatrix::LU() const {
		CSquareMatrix D{ m_Size };
    D.copyFrom( *this );

    for( auto k = 0u; k <= m_Size - 2; ++k ) {
      for( auto j = k + 1; j <= m_Size - 1; ++j) {
	      const auto x = D[ j ][ k ] / D[ k ][ k ];
        for( auto i = k; i <= m_Size - 1; ++i ) {  
          D[ j ][ i ] = D[ j ][ i ] - x * D[ k ][ i ];
        }
        D[ j ][ k ] = x;
      }
    }

    return D;
  }

}
