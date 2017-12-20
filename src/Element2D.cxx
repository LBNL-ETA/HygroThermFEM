#include "Element2D.hxx"
#include "IntegrationPoints.hxx"
#include "QuadrilateralLocal2D.hxx"

using FenestrationCommon::SquareMatrix;

namespace MoisThermFEM {

	//////////////////////////////////////////////////////////////////////////////
	//  IElementQuadrilateral2D
	//////////////////////////////////////////////////////////////////////////////

	IElementQuadrilateral2D::IElementQuadrilateral2D( const Node2D & t_Node1, const Node2D & t_Node2,
		const Node2D & t_Node3, const Node2D & t_Node4 ) : m_Element( t_Node1, t_Node2, t_Node3, t_Node4 ),
			  m_ElementNodes( t_Node1, t_Node2, t_Node3, t_Node4 ) {

	}

	std::vector<size_t>	IElementQuadrilateral2D::nodeIndexes() const {
		return m_ElementNodes.getNodeIndexes();
	}

	//////////////////////////////////////////////////////////////////////////////
	//  IQLEMatrix2D
	//////////////////////////////////////////////////////////////////////////////
	IQLEMatrix2D::IQLEMatrix2D( const double t_Value ) : 
		m_Matrix( numOfQadrilateralNodes ), m_Value( t_Value ) {
		
	}

	SquareMatrix< double >	IQLEMatrix2D::getMatrix() const	{
		return m_Matrix;
	}

	void IQLEMatrix2D::calculate() {
		const auto count = IntegrationPoints2D::Instance().count2D();

		for( unsigned i = 0; i < count; ++i ) {
			auto aMatInIntPt = calculateMatrixInIntegrationPoint( i );
			for( unsigned j = 0; j < aMatInIntPt.size(); ++j ) {
				for( unsigned k = 0; k < aMatInIntPt.size(); ++k ) {
					m_Matrix[ j ][ k ] += aMatInIntPt[ j ][ k ];
				}
			}
		}
	}

	//////////////////////////////////////////////////////////////////////////////
	//  IQLEConductance2D
	//////////////////////////////////////////////////////////////////////////////

	IQLEConductance2D::IQLEConductance2D( const Node2D & t_Node1, const Node2D & t_Node2, 
		const Node2D & t_Node3, const Node2D & t_Node4, const double t_Value) : 
		IElementQuadrilateral2D( t_Node1, t_Node2, t_Node3, t_Node4 ), IQLEMatrix2D( t_Value ) {
		
	}

	SquareMatrix< double > IQLEConductance2D::calculateMatrixInIntegrationPoint(
		const size_t t_IntegrationPointIndex) const {
		
		SquareMatrix< double > aMatrix( numOfQadrilateralNodes );

		auto DPsiDx = m_Element.DPsiDx( t_IntegrationPointIndex );
		auto DPsiDy = m_Element.DPsiDy( t_IntegrationPointIndex );
		const auto det = m_Element.det( t_IntegrationPointIndex );

		for( size_t i = 0; i < aMatrix.size(); ++i ) {
			for( size_t j = 0; j < aMatrix.size(); ++j ) {
				aMatrix[ i ][ j ] =	( DPsiDx[ i ] * DPsiDx[ j ] + DPsiDy[ i ] * DPsiDy[ j ] ) * det * m_Value;
			}
		}

		return aMatrix;
	}

	//////////////////////////////////////////////////////////////////////////////
	//  IQLECapacitance2D
	//////////////////////////////////////////////////////////////////////////////

	IQLECapacitance2D::IQLECapacitance2D( const Node2D & t_Node1, const Node2D & t_Node2, 
		const Node2D & t_Node3, const Node2D & t_Node4, const double t_Value ) :
		IElementQuadrilateral2D( t_Node1, t_Node2, t_Node3, t_Node4 ), IQLEMatrix2D( t_Value ) {
	}

	SquareMatrix< double > IQLECapacitance2D::calculateMatrixInIntegrationPoint(
		const size_t t_IntegrationPointIndex ) const {
		SquareMatrix< double > aMatrix( numOfQadrilateralNodes );

		auto & aElement = QuadrilateralLinearLocal2D::Instance();

		auto psi = aElement.VPsi( t_IntegrationPointIndex );
		const auto det = m_Element.det( t_IntegrationPointIndex );

		for( size_t i = 0; i < aMatrix.size(); ++i ) {
			for( size_t j = 0; j < aMatrix.size(); ++j ) {
				aMatrix[ i ][ j ] = psi[ i ] * psi[ j ] * det * m_Value;
			}
		}

		return aMatrix;
	}

	//////////////////////////////////////////////////////////////////////////////
	//  ElementThermalLinear2D
	//////////////////////////////////////////////////////////////////////////////

	ElementThermalLinear2D::ElementThermalLinear2D( const Node2D & t_Node1, const Node2D & t_Node2,
	                                  const Node2D & t_Node3, const Node2D & t_Node4,
	                                  const double t_Cond, const double t_Rho, const double t_Cp ) :
		m_Conductance( t_Node1, t_Node2, t_Node3, t_Node4, t_Cond ),
		m_Capacitance( t_Node1, t_Node2, t_Node3, t_Node4, t_Cp * t_Rho )	{
		m_Capacitance.calculate();
		m_Conductance.calculate();
	}

	ElementThermalLinear2D::ElementThermalLinear2D( const ElementThermalLinear2D & t_Element ) : 
		m_Conductance( t_Element.m_Conductance ), m_Capacitance( t_Element.m_Capacitance ) {
		
	}

	std::vector< size_t > ElementThermalLinear2D::nodeIndexes() const {
		return m_Conductance.nodeIndexes();
	}

	SquareMatrix< double > ElementThermalLinear2D::conductivity() const {
		return m_Conductance.getMatrix();
	}

	SquareMatrix< double > ElementThermalLinear2D::rhoCp() const {
		return m_Capacitance.getMatrix();
	}

}
