#include <math.h>
#include <assert.h>

#include "IBCLine2D.hxx"
#include "Node2D.hxx"
#include "IntegrationPoints.hxx"
#include "LineLocal1D.hxx"

using FenestrationCommon::SquareMatrix;

namespace MoisThermFEM {

	const std::size_t numOfNodes = 2;

	IBCLinear2D::IBCLinear2D( const Node2D & t_Node1, const Node2D & t_Node2,
														const bool t_Linear )
			: m_Nodes( t_Node1, t_Node2 ), m_Linear( t_Linear ), m_PsiPsiMatrix( numOfNodes ),
				m_PsiVector( numOfNodes, 0 ) {

		m_Determinant =
				0.5 * sqrt( pow( t_Node1.X() - t_Node2.X(), 2 ) + pow( t_Node1.Y() - t_Node2.Y(), 2 ) );

		// Create psi matrix and vector
		for ( std::size_t i = 0; i < numOfIntegrationPoints(); ++i ) {
			for ( std::size_t j = 0; j < numOfNodes; ++j ) {
				for ( std::size_t k = 0; k < numOfNodes; ++k ) {
					m_PsiPsiMatrix[ j ][ k ] += m_Determinant * psi( i, j ) * psi( i, k );
				}
				m_PsiVector[ j ] += m_Determinant * psi( i, j );
			}
		}

	}

	std::vector< size_t > IBCLinear2D::getNodeIndexes() const {
		return m_Nodes.getNodeIndexes();
	}

	size_t IBCLinear2D::numOfIntegrationPoints() {
		return IntegrationPoints2D::Instance().count1D();
	}

	double IBCLinear2D::psi( const size_t IntegrationPointIndex, const size_t Index ) {
		return LineLinearLocal1D::Instance().Psi( IntegrationPointIndex, Index );
	}

	bool IBCLinear2D::isLinear() const {
		return m_Linear;
	}

	FenestrationCommon::SquareMatrix< double > IBCLinear2D::D_HMatrix() const {
		//TODO: Return zero matrix for now
		return SquareMatrix< double >( 2 );
	}

}