#include "Element2D.hxx"
#include "IntegrationPoints.hxx"
#include "QuadrilateralLocal2D.hxx"
#include "Material.hxx"
#include "FEMunique.hxx"

using FenestrationCommon::SquareMatrix;

namespace MoisThermFEM {

	//////////////////////////////////////////////////////////////////////////////
	///  IElementQuadrilateral2D
	//////////////////////////////////////////////////////////////////////////////

	IElementQuadrilateral2D::IElementQuadrilateral2D( const Node2D & t_Node1, const Node2D & t_Node2,
																										const Node2D & t_Node3, const Node2D & t_Node4 )
			: m_Element( t_Node1, t_Node2, t_Node3, t_Node4 ),
				m_ElementNodes( t_Node1, t_Node2, t_Node3, t_Node4 ) {

	}

	std::vector< size_t > IElementQuadrilateral2D::nodeIndexes() const {
		return m_ElementNodes.getNodeIndexes();
	}

	//////////////////////////////////////////////////////////////////////////////
	///  IQLEMatrix2D
	//////////////////////////////////////////////////////////////////////////////
	IQLEMatrix2D::IQLEMatrix2D( const double t_Value ) :
			m_Matrix( numOfQuadrilateralNodes ), m_Value( t_Value ) {

	}

	SquareMatrix< double > IQLEMatrix2D::getMatrix() const {
		return m_Matrix;
	}

	void IQLEMatrix2D::integrate() {
		const auto count = IntegrationPoints2D::Instance().count2D();

		for ( unsigned i = 0; i < count; ++i ) {
			auto aMatInIntPt = calculateMatrixInIntegrationPoint( i );
			for ( unsigned j = 0; j < aMatInIntPt.size(); ++j ) {
				for ( unsigned k = 0; k < aMatInIntPt.size(); ++k ) {
					m_Matrix[ j ][ k ] += aMatInIntPt[ j ][ k ];
				}
			}
		}
	}

	//////////////////////////////////////////////////////////////////////////////
	///  QLEConductance2D
	//////////////////////////////////////////////////////////////////////////////

	QLEConductance2D::QLEConductance2D( const Node2D & t_Node1, const Node2D & t_Node2,
																			const Node2D & t_Node3, const Node2D & t_Node4,
																			const double t_Value ) :
			IElementQuadrilateral2D( t_Node1, t_Node2, t_Node3, t_Node4 ), IQLEMatrix2D( t_Value ) {

	}

	SquareMatrix< double > QLEConductance2D::calculateMatrixInIntegrationPoint(
			const size_t t_IntegrationPointIndex ) const {

		SquareMatrix< double > aMatrix( numOfQuadrilateralNodes );

		auto DPsiDx = m_Element.DPsiDx( t_IntegrationPointIndex );
		auto DPsiDy = m_Element.DPsiDy( t_IntegrationPointIndex );
		const auto det = m_Element.det( t_IntegrationPointIndex );

		for ( size_t i = 0; i < aMatrix.size(); ++i ) {
			for ( size_t j = 0; j < aMatrix.size(); ++j ) {
				aMatrix[ i ][ j ] =
						( DPsiDx[ i ] * DPsiDx[ j ] + DPsiDy[ i ] * DPsiDy[ j ] ) * det * m_Value;
			}
		}

		return aMatrix;
	}

	//////////////////////////////////////////////////////////////////////////////
	//  QLECapacitance2D
	//////////////////////////////////////////////////////////////////////////////

	QLECapacitance2D::QLECapacitance2D( const Node2D & t_Node1, const Node2D & t_Node2,
																			const Node2D & t_Node3, const Node2D & t_Node4,
																			const double t_Value ) :
			IElementQuadrilateral2D( t_Node1, t_Node2, t_Node3, t_Node4 ), IQLEMatrix2D( t_Value ) {
	}

	SquareMatrix< double > QLECapacitance2D::calculateMatrixInIntegrationPoint(
			const size_t t_IntegrationPointIndex ) const {
		SquareMatrix< double > aMatrix( numOfQuadrilateralNodes );

		auto & aElement = QuadrilateralLinearLocal2D::Instance();

		auto psi = aElement.VPsi( t_IntegrationPointIndex );
		const auto det = m_Element.det( t_IntegrationPointIndex );

		for ( size_t i = 0; i < aMatrix.size(); ++i ) {
			for ( size_t j = 0; j < aMatrix.size(); ++j ) {
				aMatrix[ i ][ j ] = psi[ i ] * psi[ j ] * det * m_Value;
			}
		}

		return aMatrix;
	}

	//////////////////////////////////////////////////////////////////////////////
	///  IElementlLinear2D
	//////////////////////////////////////////////////////////////////////////////

	IElementLinear2D::IElementLinear2D(
			const Node2D & t_Node1, const Node2D & t_Node2,
			const Node2D & t_Node3, const Node2D & t_Node4,
			std::unique_ptr< FenestrationCommon::ICurve > t_Conductance,
			std::unique_ptr< FenestrationCommon::ICurve > t_Capacitance ) :
			m_Node1( t_Node1 ), m_Node2( t_Node2 ), m_Node3( t_Node3 ), m_Node4( t_Node4 ),
			m_Conductance{ std::move( t_Conductance ) }, m_Capacitance{ std::move( t_Capacitance ) } {

	}

	std::vector< size_t > IElementLinear2D::nodeIndexes() const {
		QLEConductance2D aMatrix( m_Node1, m_Node2, m_Node3, m_Node4, m_Conductance->value() );
		return aMatrix.nodeIndexes();
	}

	SquareMatrix< double > IElementLinear2D::conductanceMatrix() const {
		QLEConductance2D aMatrix( m_Node1, m_Node2, m_Node3, m_Node4, m_Conductance->value() );
		aMatrix.integrate();
		return aMatrix.getMatrix();
	}

	SquareMatrix< double > IElementLinear2D::capacitanceMatrix() const {
		QLECapacitance2D aMatrix( m_Node1, m_Node2, m_Node3, m_Node4, m_Capacitance->value() );
		aMatrix.integrate();
		return aMatrix.getMatrix();
	}

	//////////////////////////////////////////////////////////////////////////////
	///  ElementThermalLinear2D
	//////////////////////////////////////////////////////////////////////////////

	ElementThermalLinear2D::ElementThermalLinear2D( const Node2D & t_Node1, const Node2D & t_Node2,
																									const Node2D & t_Node3, const Node2D & t_Node4,
																									const Material & mat ) :
			IElementLinear2D( t_Node1, t_Node2, t_Node3, t_Node4,
												fem::make_unique< FenestrationCommon::Constant >(
														mat.thermalConductivity() ),
												fem::make_unique< FenestrationCommon::Constant >(
														mat.heatCapacity() * mat.density() ) ) {

	}

}
