#include <cassert>
#include <iostream>

#include "Element2D.hxx"
#include "IntegrationPoints.hxx"
#include "QuadrilateralLocal2D.hxx"

using FenestrationCommon::SquareMatrix;

namespace MoisThermFEM {

	using pValue = std::shared_ptr< MoisThermFEM::IValue >;

	//////////////////////////////////////////////////////////////////////////////
	///  IQLEMatrix2D
	//////////////////////////////////////////////////////////////////////////////
	IQLEMatrix2D::IQLEMatrix2D( const QuadrilateralLinearGlobal2D & t_Element,
															const double t_Value1, const double t_Value2, const double t_Value3,
															const double t_Value4 ) :
			m_Matrix{ numOfQuadrilateralNodes }, m_Values{ t_Value1, t_Value2, t_Value3, t_Value4 },
			m_Global2D{ t_Element },
			m_IntegrationMatrix{ numOfQuadrilateralNodes, SquareMatrix< double >{ numOfQuadrilateralNodes } } {

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

	QLEConductance2D::QLEConductance2D( const QuadrilateralLinearGlobal2D & t_Element,
																			const double t_Value1, const double t_Value2,
																			const double t_Value3, const double t_Value4 ) :
			IQLEMatrix2D{ t_Element, t_Value1, t_Value2, t_Value3, t_Value4 } {

		const auto numOfIntegrationPoints = IntegrationPoints2D::Instance().count2D();

		for ( std::size_t integrationPoint = 0; integrationPoint < numOfIntegrationPoints;
					++integrationPoint ) {
			auto DPsiDx = m_Global2D.DPsiDx( integrationPoint );
			auto DPsiDy = m_Global2D.DPsiDy( integrationPoint );
			const auto det = m_Global2D.det( integrationPoint );

			auto & DPsiDxDyMatrix = m_IntegrationMatrix[ integrationPoint ];
			for ( auto i = 0u; i < DPsiDxDyMatrix.size(); ++i ) {
				for ( auto j = 0u; j < DPsiDxDyMatrix.size(); ++j ) {
					DPsiDxDyMatrix[ i ][ j ] =
							( DPsiDx[ i ] * DPsiDx[ j ] + DPsiDy[ i ] * DPsiDy[ j ] ) * det;
				}
			}

		}
	}

	SquareMatrix< double > QLEConductance2D::calculateMatrixInIntegrationPoint(
			const size_t t_IntegrationPointIndex ) const {

		SquareMatrix< double > aMatrix( numOfQuadrilateralNodes );

		auto & DPsiDxDy = m_IntegrationMatrix[ t_IntegrationPointIndex ];

		for ( size_t i = 0; i < aMatrix.size(); ++i ) {
			for ( size_t j = 0; j < aMatrix.size(); ++j ) {
				aMatrix[ i ][ j ] = DPsiDxDy[ i ][ j ] * 0.5 * ( m_Values[ i ] + m_Values[ j ] );
			}
		}

		return aMatrix;
	}

//////////////////////////////////////////////////////////////////////////////
///  QLECapacitance2D
//////////////////////////////////////////////////////////////////////////////

	QLECapacitance2D::QLECapacitance2D( const QuadrilateralLinearGlobal2D & t_Element,
																			const double t_Value1, const double t_Value2,
																			const double t_Value3,
																			const double t_Value4 ) :
			IQLEMatrix2D{ t_Element, t_Value1, t_Value2, t_Value3, t_Value4 } {

		const auto numOfIntegrationPoints = IntegrationPoints2D::Instance().count2D();
		auto & aElement = QuadrilateralLinearLocal2D::Instance();

		for ( std::size_t integrationPoint = 0; integrationPoint < numOfIntegrationPoints;
					++integrationPoint ) {
			const auto & psi = aElement.VPsi( integrationPoint );
			const auto det = m_Global2D.det( integrationPoint );

			auto & psiPsiMatrix = m_IntegrationMatrix[ integrationPoint ];
			for ( auto i = 0u; i < numOfIntegrationPoints; ++i ) {
				for ( auto j = 0u; j < numOfIntegrationPoints; ++j ) {
					psiPsiMatrix[ i ][ j ] = det * psi[ i ] * psi[ j ];
				}
			}

		}

	}

	SquareMatrix< double > QLECapacitance2D::calculateMatrixInIntegrationPoint(
			const size_t t_IntegrationPointIndex ) const {

		SquareMatrix< double > aMatrix{ numOfQuadrilateralNodes };

		auto & psiPsi = m_IntegrationMatrix[ t_IntegrationPointIndex ];

		for ( size_t i = 0; i < aMatrix.size(); ++i ) {
			for ( size_t j = 0; j < numOfQuadrilateralNodes; ++j ) {
				aMatrix[ i ][ j ] = psiPsi[ i ][ j ] * 0.5 * ( m_Values[ i ] + m_Values[ j ] );
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
			m_Node{ { t_Node1, t_Node2, t_Node3, t_Node4 } },
			m_Global2D{ t_Node1, t_Node2, t_Node3, t_Node4 } {
	}

	SquareMatrix< double > IElementLinear2D::conductanceMatrix() const {
		FenestrationCommon::SquareMatrix< double > result{ numOfQuadrilateralNodes };
		for ( const auto & cond : m_Conductance ) {
			auto value1 = cond->value( m_Node[ 0 ].getState() );
			auto value2 = cond->value( m_Node[ 1 ].getState() );
			auto value3 = cond->value( m_Node[ 2 ].getState() );
			auto value4 = cond->value( m_Node[ 3 ].getState() );
			QLEConductance2D aMatrix( m_Global2D, value1, value2, value3, value4 );
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
			QLECapacitance2D aMatrix( m_Global2D, value1, value2, value3, value4 );
			aMatrix.integrate();
			result += aMatrix.getMatrix();
		}
		return result;
	}

	Node2D & IElementLinear2D::getNode( const std::size_t index ) {
		assert( index < m_Node.size() );
		return m_Node[ index ];
	}

	std::vector< std::size_t > IElementLinear2D::nodeIndexes() const {
		return m_Global2D.nodeIndexes();
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
		m_Capacitance.push_back( rhoCp );
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
				std::make_shared< MoisThermFEM::TabularFunction >( mat.sorptionCurve(),
																																	 Property::humidity ) );
		pValue sorptionDerivative{ std::make_shared< MoisThermFEM::Derivative >( sorptionCurve ) };

		m_Capacitance.push_back( sorptionDerivative );
	}

}
