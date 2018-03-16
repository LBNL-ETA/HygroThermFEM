#include <cassert>
#include <iostream>

#include "Element2D.hxx"
#include "IntegrationPoints.hxx"
#include "QuadrilateralLocal2D.hxx"
#include "Common.hxx"

using FenestrationCommon::SquareMatrix;

namespace MoisThermFEM {

	using pValue = std::shared_ptr< MoisThermFEM::IValue >;

	//////////////////////////////////////////////////////////////////////////////
	///  IQLEMatrix2D
	//////////////////////////////////////////////////////////////////////////////
	IQLEMatrix2D::IQLEMatrix2D( const QuadrilateralLinearGlobal2D & t_Element ) :
			m_Global2D{ t_Element },
			m_IntegrationMatrix{ numOfQuadrilateralNodes,
													 SquareMatrix< double >{ numOfQuadrilateralNodes } } {

	}

	FenestrationCommon::SquareMatrix< double >
	IQLEMatrix2D::integrate( const std::vector< double > & t_Values ) const {
		const auto count = IntegrationPoints2D::Instance().count2D();

		FenestrationCommon::SquareMatrix< double > aMatrix{ numOfQuadrilateralNodes };

		for ( auto i = 0u; i < count; ++i ) {
			calculateMatrixInIntegrationPoint( t_Values, i, aMatrix );
		}

		return aMatrix;
	}

	void IQLEMatrix2D::calculateMatrixInIntegrationPoint(
			const std::vector< double > & t_Values, const std::size_t t_IntegrationPointIndex,
			FenestrationCommon::SquareMatrix< double > & t_Matrix ) const {

		assert( t_Values.size() == 4 );

		auto & intPointMatrix = m_IntegrationMatrix[ t_IntegrationPointIndex ];

		for ( size_t i = 0; i < t_Matrix.size(); ++i ) {
			for ( size_t j = 0; j < t_Matrix.size(); ++j ) {
				t_Matrix[ i ][ j ] += intPointMatrix[ i ][ j ] * 0.5 * ( t_Values[ i ] + t_Values[ j ] );
			}
		}
	}

//////////////////////////////////////////////////////////////////////////////
///  QLEConductance2D
//////////////////////////////////////////////////////////////////////////////

	QLEConductance2D::QLEConductance2D( const QuadrilateralLinearGlobal2D & t_Element ) :
			IQLEMatrix2D{ t_Element } {

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

//////////////////////////////////////////////////////////////////////////////
///  QLECapacitance2D
//////////////////////////////////////////////////////////////////////////////

	QLECapacitance2D::QLECapacitance2D( const QuadrilateralLinearGlobal2D & t_Element ) :
			IQLEMatrix2D{ t_Element } {

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

//////////////////////////////////////////////////////////////////////////////
///  IElementLinear2D
//////////////////////////////////////////////////////////////////////////////

	IElementLinear2D::IElementLinear2D(
			const Node2D & t_Node1, const Node2D & t_Node2,
			const Node2D & t_Node3, const Node2D & t_Node4 ) :
			m_Node{ { t_Node1, t_Node2, t_Node3, t_Node4 } },
			m_Global2D{ t_Node1, t_Node2, t_Node3, t_Node4 },
			m_QLECapacitance2D{ m_Global2D },
			m_QLEConductance2D{ m_Global2D } {

	}

	SquareMatrix< double > IElementLinear2D::conductanceMatrix() const {
		FenestrationCommon::SquareMatrix< double > result{ numOfQuadrilateralNodes };
		const auto numOfIntegrationPoints = IntegrationPoints2D::Instance().count2D();
		for ( const auto & cond : m_Conductance ) {
			std::vector< double > values( numOfIntegrationPoints );
			for ( auto i = 0u; i < numOfIntegrationPoints; ++i ) {
				values[ i ] = cond->value( m_Node[ i ].getState() );
			}
			result += m_QLEConductance2D.integrate( values );
		}

		return result;
	}

	SquareMatrix< double > IElementLinear2D::capacitanceMatrix() const {
		FenestrationCommon::SquareMatrix< double > result{ numOfQuadrilateralNodes };
		const auto numOfIntegrationPoints = IntegrationPoints2D::Instance().count2D();
		for ( const auto & cap : m_Capacitance ) {
			std::vector< double > values( numOfIntegrationPoints );
			for ( auto i = 0u; i < numOfIntegrationPoints; ++i ) {
				values[ i ] = cap->value( m_Node[ i ].getState() );
			}
			result += m_QLECapacitance2D.integrate( values );
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

		/// Calculate air and water content
		pValue waterContent(
				std::make_shared< MoisThermFEM::TabularFunction >( mat.sorptionCurve(),
																													 Property::humidity ) );

		/// Calls sorption curve at 100% humidity to get maximum water content
		auto maxWaterContent = waterContent->value( State( 0, 1, 0 ) );

		pValue waterFill = mat.porosity() / maxWaterContent * waterContent;

		pValue airFill = mat.porosity() - waterFill;

		auto waterCapacitance = waterFill * ( Constants::Density_Water * Constants::Cp_Water );
		auto airCapacitance = airFill * ( Constants::Density_Air * Constants::Cp_Air );
		auto dryCapacitance =
				( 1 - mat.porosity() ) * ( mat.density() * mat.heatCapacity() );

		auto capacitance = waterCapacitance + airCapacitance;
		capacitance = capacitance	+ dryCapacitance;

		m_Capacitance.push_back( capacitance );

		auto waterConductance = waterFill * Constants::K_Water;
		auto airConductance = airFill * Constants::K_Air;
		auto dryConductance = ( 1 - mat.porosity() ) * mat.thermalConductivity();

		pValue thermalConductivity(
				std::make_shared< MoisThermFEM::Constant >( mat.thermalConductivity() ) );

		auto conductance = waterConductance + airConductance;
		conductance = conductance + dryConductance;

		m_Conductance.push_back( conductance );

	}

	//////////////////////////////////////////////////////////////////////////////
	///  ElementMoistureLinear2D
	//////////////////////////////////////////////////////////////////////////////

	ElementMoistureLinear2D::ElementMoistureLinear2D( const Node2D & t_Node1, const Node2D & t_Node2,
																										const Node2D & t_Node3, const Node2D & t_Node4,
																										const Material & mat ) :
			IElementLinear2D( t_Node1, t_Node2, t_Node3, t_Node4 ) {

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
