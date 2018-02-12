#pragma once

#include <vector>
#include <functional>
#include <stdexcept>
#include <algorithm>

namespace FenestrationCommon {

	template < class T >
	class Vector : public std::vector< T > {
	public:
		Vector( const std::size_t size, const T & value ) : std::vector< T >( size, value ) {}

		Vector( const std::vector< T > & t_vector ) : std::vector< T >( t_vector ) {}

		Vector( const std::initializer_list< T > & t_vector ) : std::vector< T >( t_vector ) {}

		friend Vector< T > operator+( const Vector< T > & lhs, const Vector< T > & rhs ) {
			if( rhs.size() != lhs.size() ) {
				throw std::runtime_error( "Vectors must be identical in size." );
			}

			Vector< T > result( rhs.size(), 0.0 );
			std::transform( rhs.begin(), rhs.end(), lhs.begin(), result.begin(),
											std::plus< T >() );

			return result;
		}

		friend Vector< T > operator-( const Vector< T > & lhs, const Vector< T > & rhs ) {
			if( rhs.size() != lhs.size() ) {
				throw std::runtime_error( "Vectors must be identical in size." );
			}

			Vector< T > result( rhs.size(), 0 );
			std::transform( lhs.begin(), lhs.end(), rhs.begin(), result.begin(),
											std::minus< T >() );

			return result;
		}

		friend Vector< T > operator*( const Vector< T > & lhs, const Vector< T > & rhs ) {
			if( rhs.size() != lhs.size() ) {
				throw std::runtime_error( "Vectors must be identical in size." );
			}

			Vector< T > result( rhs.size(), 0 );
			std::transform( lhs.begin(), lhs.end(), rhs.begin(), result.begin(),
											std::multiplies< T >() );

			return result;
		}

		friend Vector< T > operator*( const Vector< T > & lhs, const T value ) {

			Vector< T > result( lhs.size(), 0 );
			std::transform( lhs.begin(), lhs.end(), result.begin(),
											std::bind1st( std::multiplies< T >(), value ) );

			return result;
		}

		friend Vector< T > operator*( const T value, const Vector< T > & lhs ) {
			return lhs * value;
		}

	};

}