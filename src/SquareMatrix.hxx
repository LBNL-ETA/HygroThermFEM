#pragma once

#include <vector>
#include <stdexcept>

#include "Vector.hxx"

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

		explicit SquareMatrix( const std::initializer_list< std::vector< T > > & t_input ) : m_Matrix(
				t_input ), m_size( m_Matrix.size() ) {
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

		SquareMatrix addDiagonal( const std::vector< T > & t_Vector ) {
			if( m_size != t_Vector.size() ) {
				std::runtime_error( "Matrix and vector have different sizes." );
			}

			SquareMatrix aMatrix { m_size };
			for ( auto i = 0u; i < m_size; ++i ) {
				for ( auto j = 0u; j < m_size; ++j ) {
					aMatrix[ i ][ j ] = m_Matrix[ i ][ j ];
				}
				aMatrix[ i ][ i ] += t_Vector[ i ];
			}

			return aMatrix;
		}

		SquareMatrix< T > &
		operator+( const SquareMatrix< T > & rhs ) {
			if( rhs.m_size != m_size ) {
				throw std::runtime_error( "Matrices must be identical in size." );
			}

			for ( unsigned i = 0; i < m_size; ++i ) {
				for ( unsigned j = 0; j < m_size; ++j ) {
					m_Matrix[ i ][ j ] = m_Matrix[ i ][ j ] + rhs.m_Matrix[ i ][ j ];
				}
			}

			return *this;
		}

		SquareMatrix< T > &
		operator-( const SquareMatrix< T > & rhs ) {
			if( rhs.m_size != m_size ) {
				throw std::runtime_error( "Matrices must be identical in size." );
			}

			for ( unsigned i = 0; i < m_size; ++i ) {
				for ( unsigned j = 0; j < m_size; ++j ) {
					m_Matrix[ i ][ j ] = m_Matrix[ i ][ j ] - rhs.m_Matrix[ i ][ j ];
				}
			}

			return *this;
		}

		SquareMatrix< T > operator*( const SquareMatrix< T > & rhs ) {
			if( m_size != rhs.size() ) {
				throw std::runtime_error( "Matrices must be identical in size." );
			}

			SquareMatrix< T > temp( m_size );

			for ( auto i = 0u; i < m_size; ++i ) {
				for ( auto k = 0u; k < m_size; ++k ) {
					for ( auto j = 0u; j < m_size; ++j ) {
						temp[ i ][ j ] = temp[ i ][ j ] + m_Matrix[ i ][ k ] * rhs.m_Matrix[ k ][ j ];
					}
				}
			}

			return temp;
		}

		friend SquareMatrix< T >
		operator*( const SquareMatrix< T > & lhs, const double t_Value ) {

			SquareMatrix aMatrix { lhs.size() };

			for ( auto i = 0u; i < lhs.size(); ++i ) {
				for ( auto k = 0u; k < lhs.size(); ++k ) {
					aMatrix[ i ][ k ] += t_Value * lhs.m_Matrix[ i ][ k ];
				}
			}

			return aMatrix;
		}

		friend SquareMatrix< T >
		operator*( const double t_Value, const SquareMatrix< T > & lhs ) {
			return lhs * t_Value;
		}

		friend FenestrationCommon::Vector< T >
		operator*( const FenestrationCommon::Vector< T > & t_vector, const SquareMatrix & t_matrix ) {
			if( t_vector.size() != t_matrix.size() ) {
				throw std::runtime_error( "Vector and matrix have incompatible sizes." );
			}

			std::vector< T > multV( t_vector.size() );
			for ( auto i = 0u; i < t_vector.size(); ++i ) {
				for ( auto j = 0u; j < t_vector.size(); ++j ) {
					multV[ i ] += t_matrix.m_Matrix[ j ][ i ] * t_vector[ j ];
				}
			}

			return multV;
		}

		friend FenestrationCommon::Vector< T >
		operator*( const SquareMatrix & t_matrix, const FenestrationCommon::Vector< T > & t_vector ) {
			if( t_vector.size() != t_matrix.size() ) {
				throw std::runtime_error( "Vector and matrix have incompatible sizes." );
			}

			std::vector< T > multV( t_vector.size() );
			for ( auto i = 0u; i < t_vector.size(); ++i ) {
				for ( auto j = 0u; j < t_vector.size(); ++j ) {
					multV[ i ] += t_matrix.m_Matrix[ i ][ j ] * t_vector[ j ];
				}
			}

			return multV;
		}

	protected:
		std::vector< std::vector< T > > m_Matrix;
		std::size_t m_size;
	};

}