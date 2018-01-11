#include <math.h>
#include <assert.h>

#include "Line2D.hxx"
#include "Node2D.hxx"
#include "IntegrationPoints.hxx"
#include "LineLocal1D.hxx"

using FenestrationCommon::SquareMatrix;

namespace MoisThermFEM {

	ILineLinear2D::ILineLinear2D( const Node2D & t_Node1, const Node2D & t_Node2 )
			: m_Nodes( t_Node1, t_Node2 ), m_Rvector( 2, 0 ), m_matrixA( 2 ) {

		m_Determinant =
				0.5 * sqrt( pow( t_Node1.X() - t_Node2.X(), 2 ) + pow( t_Node1.Y() - t_Node2.Y(), 2 ) );

	}

	std::vector< size_t > ILineLinear2D::getNodeIndexes() const {
		return m_Nodes.getNodeIndexes();
	}

	std::vector< double > ILineLinear2D::rightHandSideVector() const {
		return m_Rvector;
	}

	SquareMatrix< double > ILineLinear2D::matrixA() const {
		return m_matrixA;
	}

	size_t ILineLinear2D::numOfIntegrationPoints() {
		return IntegrationPoints2D::Instance().count1D();
	}

	double ILineLinear2D::psi( const size_t IntegrationPointIndex, const size_t Index ) {
		return LineLinearLocal1D::Instance().Psi( IntegrationPointIndex, Index );
	}

}