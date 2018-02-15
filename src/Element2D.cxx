#include "Element2D.hxx"
#include "IntegrationPoints.hxx"
#include "QuadrilateralLocal2D.hxx"
#include "FEMunique.hxx"

using FenestrationCommon::SquareMatrix;

namespace MoisThermFEM {

	//////////////////////////////////////////////////////////////////////////////
	///  IElementQuadrilateral2D
	//////////////////////////////////////////////////////////////////////////////

	IElementQuadrilateral2D::IElementQuadrilateral2D( const Node2D & t_Node1, const Node2D & t_Node2,
																										const Node2D & t_Node3, const Node2D & t_Node4 )
			: m_Element( t_Node1, t_Node2, t_Node3, t_Node4 ) {

	}

	std::vector< size_t > IElementQuadrilateral2D::nodeIndexes() const {
		return m_Element.nodeIndexes();
	}

	//////////////////////////////////////////////////////////////////////////////
	///  IQLEMatrix2D
	//////////////////////////////////////////////////////////////////////////////
	IQLEMatrix2D::IQLEMatrix2D( const double t_Value1, const double t_Value2, const double t_Value3,
															const double t_Value4 ) :
			m_Matrix( numOfQuadrilateralNodes ), m_Values{ t_Value1, t_Value2, t_Value3, t_Value4 } {

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
																			const double t_Value1, const double t_Value2,
																			const double t_Value3,
																			const double t_Value4 ) :
			IElementQuadrilateral2D{ t_Node1, t_Node2, t_Node3, t_Node4 },
			IQLEMatrix2D{ t_Value1, t_Value2, t_Value3, t_Value4 } {

	}

	SquareMatrix< double > QLEConductance2D::calculateMatrixInIntegrationPoint(
			const size_t t_IntegrationPointIndex ) const {

		SquareMatrix< double > aMatrix( numOfQuadrilateralNodes );

		auto DPsiDx = m_Element.DPsiDx( t_IntegrationPointIndex );
		auto DPsiDy = m_Element.DPsiDy( t_IntegrationPointIndex );
		const auto det = m_Element.det( t_IntegrationPointIndex );

		for ( size_t i = 0; i < aMatrix.size(); ++i ) {
			for ( size_t j = 0; j < aMatrix.size(); ++j ) {
				aMatrix[ i ][ j ] = ( DPsiDx[ i ] * DPsiDx[ j ] + DPsiDy[ i ] * DPsiDy[ j ] ) * det *
														( m_Values[ i ] + m_Values[ j ] ) / 2;
			}
		}

		return aMatrix;
	}

	//////////////////////////////////////////////////////////////////////////////
	//  QLECapacitance2D
	//////////////////////////////////////////////////////////////////////////////

	QLECapacitance2D::QLECapacitance2D( const Node2D & t_Node1, const Node2D & t_Node2,
																			const Node2D & t_Node3, const Node2D & t_Node4,
																			const double t_Value1, const double t_Value2,
																			const double t_Value3,
																			const double t_Value4 ) :
			IElementQuadrilateral2D{ t_Node1, t_Node2, t_Node3, t_Node4 },
			IQLEMatrix2D{ t_Value1, t_Value2, t_Value3, t_Value4 } {
	}

	SquareMatrix< double > QLECapacitance2D::calculateMatrixInIntegrationPoint(
			const size_t t_IntegrationPointIndex ) const {
		SquareMatrix< double > aMatrix( numOfQuadrilateralNodes );

		auto & aElement = QuadrilateralLinearLocal2D::Instance();

		auto psi = aElement.VPsi( t_IntegrationPointIndex );
		const auto det = m_Element.det( t_IntegrationPointIndex );

		for ( size_t i = 0; i < aMatrix.size(); ++i ) {
			for ( size_t j = 0; j < aMatrix.size(); ++j ) {
				aMatrix[ i ][ j ] = psi[ i ] * psi[ j ] * det * ( m_Values[ i ] + m_Values[ j ] ) / 2;
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
			const Property t_property ) : IElementQuadrilateral2D( t_Node1, t_Node2, t_Node3, t_Node4 ),
																		m_Node1( t_Node1 ), m_Node2( t_Node2 ), m_Node3( t_Node3 ),
																		m_Node4( t_Node4 ),
																		m_Property( t_property ) {

	}

	SquareMatrix< double > IElementLinear2D::conductanceMatrix() const {
		FenestrationCommon::SquareMatrix< double > result{ numOfQuadrilateralNodes };
		for ( const std::unique_ptr< MoisThermFEM::IFunction > & cond : m_Conductance ) {
			auto value1 = cond->value( m_Node1.getState() );
			auto value2 = cond->value( m_Node2.getState() );
			auto value3 = cond->value( m_Node3.getState() );
			auto value4 = cond->value( m_Node4.getState() );
			QLEConductance2D aMatrix( m_Node1, m_Node2, m_Node3, m_Node4, value1, value2, value3,
																value4 );
			aMatrix.integrate();
			result += aMatrix.getMatrix();
		}

		return result;
	}

	SquareMatrix< double > IElementLinear2D::capacitanceMatrix() const {
		FenestrationCommon::SquareMatrix< double > result{ numOfQuadrilateralNodes };
		for ( const std::unique_ptr< MoisThermFEM::IFunction > & cap : m_Capacitance ) {
			auto value1 = cap->value( m_Node1.getState() );
			auto value2 = cap->value( m_Node2.getState() );
			auto value3 = cap->value( m_Node3.getState() );
			auto value4 = cap->value( m_Node4.getState() );
			QLECapacitance2D aMatrix( m_Node1, m_Node2, m_Node3, m_Node4, value1, value2, value3,
																value4 );
			aMatrix.integrate();
			result += aMatrix.getMatrix();
		}
		return result;
	}

	//////////////////////////////////////////////////////////////////////////////
	///  ElementThermalLinear2D
	//////////////////////////////////////////////////////////////////////////////

	ElementThermalLinear2D::ElementThermalLinear2D( const Node2D & t_Node1, const Node2D & t_Node2,
																									const Node2D & t_Node3, const Node2D & t_Node4,
																									const Material & mat ) :
			IElementLinear2D( t_Node1, t_Node2, t_Node3, t_Node4, Property::temperature ) {
		/// Note that this works for non-porous material with constant properties.
		m_Conductance.push_back(
				fem::make_unique< MoisThermFEM::Constant >( mat.thermalConductivity(),
																										Property::temperature ) );
		m_Capacitance.push_back(
				fem::make_unique< MoisThermFEM::Constant >( mat.heatCapacity() * mat.density(),
																										Property::temperature ) );
	}

}
