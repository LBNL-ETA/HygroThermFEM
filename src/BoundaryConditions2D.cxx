#include "BoundaryConditions2D.hxx"
#include "Line2D.hxx"
#include "NodePool.hxx"
#include "SquareMatrix.hxx"

using namespace FenestrationCommon;

namespace MoisThermFEM {

	BoundaryConditions2D::BoundaryConditions2D(
			std::vector< std::unique_ptr< ILineLinear2D > > const & t_BCs ) : m_MatrixA(
			NodePool::Instance().maxIndex() ), m_vectorR( NodePool::Instance().maxIndex() ) {

		// Create full size matrices
		for( const auto & aBc : t_BCs ) {
			auto indexes = aBc->getNodeIndexes();
			auto matA = aBc->matrixA();
			auto vecR = aBc->rightHandSideVector();
			for( size_t i = 0; i < 2; ++i ) {
				for( size_t j = 0; j < 2; ++j ) {
					m_MatrixA[ indexes[ i ] - 1 ][ indexes[ j ] - 1 ] += matA[ i ][ j ];
				}
				m_vectorR[ indexes[ i ] - 1 ] += vecR[ i ];
			}
		}
	}

	BoundaryConditions2D::BoundaryConditions2D( const BoundaryConditions2D & other ) : m_MatrixA(
			NodePool::Instance().maxIndex() ) {
		*this = other;
	}

	BoundaryConditions2D & BoundaryConditions2D::operator=( const BoundaryConditions2D & other ) {
		m_vectorR = other.m_vectorR;
		m_MatrixA = other.m_MatrixA;

		return *this;
	}

	CSquareMatrix BoundaryConditions2D::matrixA() const {
		return m_MatrixA;
	}

	std::vector< double > BoundaryConditions2D::vectorR() const {
		return m_vectorR;
	}

}