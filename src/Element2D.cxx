#include <cassert>

#include "Element2D.hxx"
#include "IntegrationPoints.hxx"
#include "QuadrilateralLocal2D.hxx"
#include "FEMunique.hxx"

using FenestrationCommon::SquareMatrix;

namespace MoisThermFEM {

	using pValue = std::shared_ptr< MoisThermFEM::IValue >;

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
	///  QLECapacitance2D
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
	///  IElementLinear2D
	//////////////////////////////////////////////////////////////////////////////

	IElementLinear2D::IElementLinear2D(
			const Node2D & t_Node1, const Node2D & t_Node2,
			const Node2D & t_Node3, const Node2D & t_Node4 ) :
			IElementQuadrilateral2D( t_Node1, t_Node2, t_Node3, t_Node4 ),
			m_Node( { t_Node1, t_Node2, t_Node3, t_Node4 } ) {
	}

	SquareMatrix< double > IElementLinear2D::conductanceMatrix() const {
		FenestrationCommon::SquareMatrix< double > result{ numOfQuadrilateralNodes };
		for ( const auto & cond : m_Conductance ) {
			auto value1 = cond->value( m_Node[ 0 ].getState() );
			auto value2 = cond->value( m_Node[ 1 ].getState() );
			auto value3 = cond->value( m_Node[ 2 ].getState() );
			auto value4 = cond->value( m_Node[ 3 ].getState() );
			QLEConductance2D aMatrix( m_Node[ 0 ], m_Node[ 1 ], m_Node[ 2 ], m_Node[ 3 ], value1, value2,
																value3,
																value4 );
			aMatrix.integrate();
			result += aMatrix.getMatrix();
		}

		return result;
	}

	SquareMatrix< double > IElementLinear2D::capacitanceMatrix() const {
		FenestrationCommon::SquareMatrix< double > result{ numOfQuadrilateralNodes };
		for ( const auto & cap : m_Capacitance ) {
			auto value1 = cap->value( m_Node[ 0 ].getState() );
			auto value2 = cap->value( m_Node[ 1 ].getState() );
			auto value3 = cap->value( m_Node[ 2 ].getState() );
			auto value4 = cap->value( m_Node[ 3 ].getState() );
			QLECapacitance2D aMatrix( m_Node[ 0 ], m_Node[ 1 ], m_Node[ 2 ], m_Node[ 3 ], value1, value2,
																value3, value4 );
			aMatrix.integrate();
			result += aMatrix.getMatrix();
		}
		return result;
	}

	Node2D & IElementLinear2D::getNode( const std::size_t index ) {
		assert( index < m_Node.size() );
		return m_Node[ index ];
	}

	//////////////////////////////////////////////////////////////////////////////
	///  ElementThermalLinear2D
	//////////////////////////////////////////////////////////////////////////////

	ElementThermalLinear2D::ElementThermalLinear2D( const Node2D & t_Node1, const Node2D & t_Node2,
																									const Node2D & t_Node3, const Node2D & t_Node4,
																									const Material & mat ) :
			IElementLinear2D( t_Node1, t_Node2, t_Node3, t_Node4 ) {
		pValue thermalConductivity(
				std::make_shared< MoisThermFEM::Constant >( mat.thermalConductivity() ) );
		m_Conductance.push_back( thermalConductivity );
		pValue rhoCp(
				std::make_shared< MoisThermFEM::Constant >( mat.heatCapacity() * mat.density() ) );
		m_Capacitance.push_back(
				fem::make_unique< MoisThermFEM::Constant >( mat.heatCapacity() * mat.density() ) );
	}

	//////////////////////////////////////////////////////////////////////////////
	///  ElementMoistureLinear2D
	//////////////////////////////////////////////////////////////////////////////

	ElementMoistureLinear2D::ElementMoistureLinear2D( const Node2D & t_Node1, const Node2D & t_Node2,
																										const Node2D & t_Node3, const Node2D & t_Node4,
																										const Material & mat ) :
			IElementLinear2D( t_Node1, t_Node2, t_Node3, t_Node4 ) {

		//////////////////////////////////////////////////////////////////////////////
		/// Creating conductance function for vapor
		//////////////////////////////////////////////////////////////////////////////

		//pValue waterContent(
		//		std::make_shared< MoisThermFEM::TabularFunction >( mat.sorptionCurve(),
		//																											 Property::humidity ) );

		/// Calls sorption curve at 100% humidity to get maximum water content
		///auto maxWaterContent = waterContent->value( State( 0, 1, 0 ) );

		//pValue waterFill = mat.porosity() / maxWaterContent * waterContent;

		///pValue airFill = mat.porosity() - waterFill;

		pValue saturationFunction(
				std::make_shared< MoisThermFEM::SaturationFunction >( Property::temperature ) );

		m_Conductance.push_back( 2.5E-5 / mat.diffusionResistanceFactor() * saturationFunction );

		//////////////////////////////////////////////////////////////////////////////
		/// Creating conductance function for liquid
		//////////////////////////////////////////////////////////////////////////////
		pValue suctionCurve(
				std::make_shared< MoisThermFEM::SuctionFunction >( mat.liquidTransportationCurve(),
																													 Property::humidity ) );
		m_Conductance.push_back( suctionCurve );

		//////////////////////////////////////////////////////////////////////////////
		/// Creating capacitance function
		//////////////////////////////////////////////////////////////////////////////
		pValue sorptionCurve(
				std::make_shared< MoisThermFEM::FirstDerivativeFunction >( mat.sorptionCurve(),
																																	 Property::humidity ) );
		m_Capacitance.push_back( sorptionCurve );
	}
}
