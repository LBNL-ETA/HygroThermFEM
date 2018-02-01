#pragma once

#include <vector>
#include <stdexcept>

namespace FenestrationCommon {

	//////////////////////////////////////////////////////////////////////////////
	// Square matrices
	//////////////////////////////////////////////////////////////////////////////
	template < class T >
	class SquareMatrix {
	public:
		explicit SquareMatrix( const std::size_t size ) :
				m_Matrix( size, std::vector< T >( size, 0 ) ), m_size( size ) {

		}

		std::vector< std::vector< T > > getMatrix() const {
			return m_Matrix;
		}

		std::vector< T > & operator[]( const std::size_t index ) {
			return m_Matrix[ index ];
		}

		std::size_t size() const {
			return m_size;
		}

		void setZeros() {
			for ( unsigned i = 0; i < m_size; ++i ) {
				for ( unsigned j = 0; j < m_size; ++j ) {
					m_Matrix[ i ][ j ] = 0;
				}
			}
		}

		void setIdentity() {
			setZeros();
			for ( unsigned i = 0; i < m_size; ++i ) {
				m_Matrix[ i ][ i ] = 1;
			}
		}

		void setDiagonal( const std::vector< T > & t_Values ) {
			if( t_Values.size() != m_size ) {
				throw std::runtime_error( "Supplied vector size mismatch matrix size" );
			}
			setZeros();
			for ( unsigned i = 0; i < m_size; ++i ) {
				m_Matrix[ i ][ i ] = t_Values[ i ];
			}
		}

		SquareMatrix addDiagonal( const std::vector< double > & t_Vector ) {
			if( m_size != t_Vector.size() ) {
				std::runtime_error( "Matrix and vector have different sizes." );
			}

			SquareMatrix aMatrix{ m_size };
			for ( unsigned i = 0; i < m_size; ++i ) {
				for ( unsigned j = 0; j < m_size; ++j ) {
					aMatrix[ i ][ j ] = m_Matrix[ i ][ j ];
				}
				aMatrix[ i ][ i ] += t_Vector[ i ];
			}

			return aMatrix;
		}

		SquareMatrix add( const SquareMatrix & t_Matrix ) const {
			if( m_size != t_Matrix.m_size ) {
				throw std::runtime_error( "Matrices must be identical in size." );
			}

			SquareMatrix aMatrix{ m_size };
			for ( unsigned i = 0; i < m_size; ++i ) {
				for ( unsigned j = 0; j < t_Matrix.size(); ++j ) {
					aMatrix[ i ][ j ] = m_Matrix[ i ][ j ] + t_Matrix.m_Matrix[ i ][ j ];
				}
			}

			return aMatrix;
		}

		SquareMatrix sub( const SquareMatrix & t_Matrix ) const {
			if( m_size != t_Matrix.m_size ) {
				throw std::runtime_error( "Matrices must be identical in size." );
			}

			SquareMatrix aMatrix{ m_size };
			for ( size_t i = 0; i < m_size; ++i ) {
				for ( size_t j = 0; j < t_Matrix.m_size; ++j ) {
					aMatrix[ i ][ j ] = m_Matrix[ i ][ j ] - t_Matrix.m_Matrix[ i ][ j ];
				}
			}

			return aMatrix;
		}

		SquareMatrix mult( const SquareMatrix & t_Matrix ) const {
			if( m_size != t_Matrix.m_size ) {
				throw std::runtime_error( "Matrices must be identical in size." );
			}

			SquareMatrix aMatrix{ m_size };

			for ( auto i = 0u; i < m_size; ++i ) {
				for ( auto k = 0u; k < m_size; ++k ) {
					for ( auto j = 0u; j < t_Matrix.m_size; ++j ) {
						aMatrix[ i ][ j ] =
								aMatrix[ i ][ j ] + m_Matrix[ i ][ k ] * t_Matrix.m_Matrix[ k ][ j ];
					}
				}
			}

			return aMatrix;
		}

		SquareMatrix mult( const double t_Value ) const {

			SquareMatrix aMatrix{ m_size };

			for ( auto i = 0u; i < m_size; ++i ) {
				for ( auto k = 0u; k < m_size; ++k ) {
					aMatrix[ i ][ k ] += t_Value * m_Matrix[ i ][ k ];
				}
			}

			return aMatrix;
		}

		std::vector< double > multMxV( const std::vector< double > & t_Vector ) {
			if( m_size != t_Vector.size() ) {
				throw std::runtime_error(
						"Vector and matrix have different sizes. Multiplication operation cannot be performed." );
			}
			std::vector< double > multV( m_size );
			for( auto i = 0u; i <  m_size; ++i ) {
				for( auto j = 0u; j <  m_size; ++j ) {
					multV[ i ] += m_Matrix[ i ][ j ] * t_Vector[ j ];
				}
			}

			return multV;
		}

	protected:
		std::vector< std::vector< T > > m_Matrix;
		std::size_t m_size;
	};

}