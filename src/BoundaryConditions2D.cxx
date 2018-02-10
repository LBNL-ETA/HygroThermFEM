#include "BoundaryConditions2D.hxx"
#include "IBCLine2D.hxx"
#include "NodePool.hxx"
#include "SquareMatrix.hxx"

using namespace FenestrationCommon;

namespace MoisThermFEM {

	BoundaryConditions2D::BoundaryConditions2D(
			std::vector< std::reference_wrapper< IBCLinear2D > > & t_BCs ) :
			m_BCs( t_BCs ), m_Linear( true ) {
		for ( const auto & aBc : m_BCs ) {
			m_Linear = m_Linear && aBc.get().isLinear();
		}
	}

	SquareMatrix< double > BoundaryConditions2D::HMatrix() const {
		FenestrationCommon::SquareMatrix< double > result{ NodePool::Instance().maxIndex() };
		for ( const auto & aBc : m_BCs ) {
			const auto & bc = aBc.get();
			auto indexes = bc.getNodeIndexes();
			auto matH = bc.H_Matrix();
			for ( size_t i = 0; i < 2; ++i ) {
				for ( size_t j = 0; j < 2; ++j ) {
					result[ indexes[ i ] - 1 ][ indexes[ j ] - 1 ] += matH[ i ][ j ];
				}
			}
		}
		return result;
	}

	FenestrationCommon::Vector< double > BoundaryConditions2D::RVector() const {
		FenestrationCommon::Vector< double > result( NodePool::Instance().maxIndex(), 0 );
		// Create full size matrices
		for ( const auto & aBc : m_BCs ) {
			const auto & bc = aBc.get();
			auto indexes = bc.getNodeIndexes();
			auto vecR = bc.R_Vector();
			for ( size_t i = 0; i < 2; ++i ) {
				result[ indexes[ i ] - 1 ] += vecR[ i ];
			}
		}
		return result;
	}

	bool BoundaryConditions2D::isLinear() const {
		return m_Linear;
	}

	FenestrationCommon::SquareMatrix< double > BoundaryConditions2D::DHMatrix() const {
		FenestrationCommon::SquareMatrix< double > HDMatrix{ NodePool::Instance().maxIndex() };
		for ( const auto & aBc : m_BCs ) {
			const auto & bc = aBc.get();
			auto indexes = bc.getNodeIndexes();
			auto matDH = bc.D_HMatrix();
			for ( size_t i = 0; i < numOfBCNodes; ++i ) {
				for ( size_t j = 0; j < numOfBCNodes; ++j ) {
					HDMatrix[ indexes[ i ] - 1 ][ indexes[ j ] - 1 ] += matDH[ i ][ j ];
				}
			}
		}
		return HDMatrix;
	}

	void BoundaryConditions2D::updateNodeTemperatures( const std::vector< double > & temperatures ) {
		for ( auto & aBc : m_BCs ) {
			auto & bc = aBc.get();
			for( auto i = 0u; i < numOfBCNodes; ++i ) {
				auto & node = bc.getNode( i );
				auto index = node.getNodeNumber();
				node.setProperty( Property::temperature, temperatures[ index - 1 ] );
			}

		}
	}

}