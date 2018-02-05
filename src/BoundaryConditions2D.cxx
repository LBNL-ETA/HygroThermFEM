#include "BoundaryConditions2D.hxx"
#include "IBCLine2D.hxx"
#include "NodePool.hxx"
#include "SquareMatrix.hxx"

using namespace FenestrationCommon;

namespace MoisThermFEM {

	BoundaryConditions2D::BoundaryConditions2D(
			const std::vector< std::reference_wrapper< const IBCLinear2D > > & t_BCs ) : m_HMatrix(
			NodePool::Instance().maxIndex() ), m_RVector( NodePool::Instance().maxIndex(), 0 ), m_Linear(
			true ) {

		// Create full size matrices
		for ( const auto & aBc : t_BCs ) {
			const auto & bc = aBc.get();
			m_Linear = m_Linear && !bc.isLinear();
			auto indexes = bc.getNodeIndexes();
			auto matA = bc.H_Matrix();
			auto vecR = bc.R_Vector();
			for ( size_t i = 0; i < 2; ++i ) {
				for ( size_t j = 0; j < 2; ++j ) {
					m_HMatrix[ indexes[ i ] - 1 ][ indexes[ j ] - 1 ] += matA[ i ][ j ];
				}
				m_RVector[ indexes[ i ] - 1 ] += vecR[ i ];
			}
		}
	}

	BoundaryConditions2D::BoundaryConditions2D( const BoundaryConditions2D & other ) : m_HMatrix(
			NodePool::Instance().maxIndex() ), m_RVector( 0, 0 ) {
		*this = other;
	}

	BoundaryConditions2D & BoundaryConditions2D::operator=( const BoundaryConditions2D & other ) {
		m_RVector = other.m_RVector;
		m_HMatrix = other.m_HMatrix;

		return *this;
	}

	SquareMatrix< double > BoundaryConditions2D::HMatrix() const {
		return m_HMatrix;
	}

	FenestrationCommon::Vector< double > BoundaryConditions2D::RVector() const {
		return m_RVector;
	}

	bool BoundaryConditions2D::isLinear() const {
			return m_Linear;
	}

}