#ifndef FENESTRATIONSQUAREMATRIX_H
#define FENESTRATIONSQUAREMATRIX_H

#include <vector>
#include <memory>

namespace FenestrationCommon {

  class CSquareMatrix {
  public:
    explicit CSquareMatrix( size_t const aSize );
    size_t getSize() const;
    void setZeros();
    // All diagonal items are one and all non diagonal are zero
    void setIdentity();
    // set diagonal values from vector
    void setDiagonal( std::vector< double > const & t_Values );
    std::vector< double >& operator[]( size_t const index );
    CSquareMatrix add( CSquareMatrix const & t_Matrix ) const;
    CSquareMatrix addDiagonal( std::vector< double > const & t_Vector );
		CSquareMatrix sub( CSquareMatrix const & t_Matrix ) const;
    CSquareMatrix mult( CSquareMatrix const & t_Matrix ) const;
    // Matrix multiplication with vector
    std::vector< double > multMxV( std::vector< double > const & t_Vector ) const;
    // Matrix multiplication with vector
    std::vector< double > multVxM( std::vector< double > const & t_Vector ) const;
    void copyFrom( CSquareMatrix const & t_Matrix );
    // inverse matrix
    CSquareMatrix inverse();

  private:
    // LU decomposition of current matrix
    CSquareMatrix LU() const;

    size_t m_Size;
    std::vector< std::vector < double > > m_Matrix;

  };
}


#endif