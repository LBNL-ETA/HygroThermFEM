#include <cassert>
#include <iostream>

#include "Element2D.hxx"
#include "IntegrationPoints.hxx"
#include "QuadrilateralLocal2D.hxx"
#include "Common.hxx"
#include "MaterialProperties.hxx"
#include "NodePool.hxx"

using FenestrationCommon::SparceSquareMatrix;

namespace MoisThermFEM {

	using iValue = std::shared_ptr< MoisThermFEM::IValue >;

	//////////////////////////////////////////////////////////////////////////////
	///  IQLEMatrix2D
	//////////////////////////////////////////////////////////////////////////////
	IQLEMatrix2D::IQLEMatrix2D( const QuadrilateralLinearGlobal2D & t_Element ) :
			m_Global2D{ t_Element },
			m_IntegrationMatrix{ numOfQuadrilateralNodes,
								 SparceSquareMatrix< double >{ numOfQuadrilateralNodes } } {

	}

	FenestrationCommon::SparceSquareMatrix< double >
	IQLEMatrix2D::integrate( const std::vector< double > & t_Values ) const {
		const auto count = IntegrationPoints2D::Instance().count2D();

		FenestrationCommon::SparceSquareMatrix< double > aMatrix{ numOfQuadrilateralNodes };

		for ( auto i = 0u; i < count; ++i ) {
			calculateMatrixInIntegrationPoint( t_Values, i, aMatrix );
		}

		return aMatrix;
	}

	void IQLEMatrix2D::calculateMatrixInIntegrationPoint(
			const std::vector< double > & t_Values, const std::size_t t_IntegrationPointIndex,
			FenestrationCommon::SparceSquareMatrix< double > & t_Matrix ) const {

		assert( t_Values.size() == 4 );

		auto & intPointMatrix = m_IntegrationMatrix[ t_IntegrationPointIndex ];

		for ( size_t i = 0; i < t_Matrix.size(); ++i ) {
			for ( size_t j = 0; j < t_Matrix.size(); ++j ) {
				t_Matrix( i, j ) += intPointMatrix( i, j ) * 0.5 * ( t_Values[ i ] + t_Values[ j ] );
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
					DPsiDxDyMatrix( i, j ) =
							( DPsiDx[ i ] * DPsiDx[ j ] + DPsiDy[ i ] * DPsiDy[ j ] ) * det;
				}
			}

		}
	}

//////////////////////////////////////////////////////////////////////////////
///  QLEConductanceDerivative2D
//////////////////////////////////////////////////////////////////////////////

	QLEConductanceDerivative2D::QLEConductanceDerivative2D(
			const QuadrilateralLinearGlobal2D & t_Element ) : IQLEMatrix2D{ t_Element } {

	}

	void QLEConductanceDerivative2D::updateIntegrationMatrix(
			const std::vector< double > & t_Values ) {

		const auto numOfIntegrationPoints = IntegrationPoints2D::Instance().count2D();
		auto & aElement = QuadrilateralLinearLocal2D::Instance();

		assert( t_Values.size() == numOfIntegrationPoints );

		for ( std::size_t integrationPoint = 0; integrationPoint < numOfIntegrationPoints;
					++integrationPoint ) {
			const auto & psi = aElement.VPsi( integrationPoint );
			auto DPsiDx = m_Global2D.DPsiDx( integrationPoint );
			auto DPsiDy = m_Global2D.DPsiDy( integrationPoint );
			const auto det = m_Global2D.det( integrationPoint );

			auto gammaX = 0.0;
			auto gammaY = 0.0;
			for ( auto k = 0u; k < numOfIntegrationPoints; ++k ) {
				gammaX += DPsiDx[ k ] * t_Values[ k ];
				gammaY += DPsiDy[ k ] * t_Values[ k ];
			}

			auto & psiPsiMatrix = m_IntegrationMatrix[ integrationPoint ];
			for ( auto i = 0u; i < numOfIntegrationPoints; ++i ) {
				for ( auto j = 0u; j < numOfIntegrationPoints; ++j ) {
					psiPsiMatrix(i, j ) =
							det * ( DPsiDx[ i ] * psi[ j ] * gammaX + DPsiDy[ i ] * psi[ j ] * gammaY );
				}
			}

		}
	}

	void QLEConductanceDerivative2D::clearIntegrationMatrix() {
		m_IntegrationMatrix.clear();
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
					psiPsiMatrix(i, j) = det * psi[ i ] * psi[ j ];
				}
			}

		}

	}

//////////////////////////////////////////////////////////////////////////////
///  DerivativeFunction
//////////////////////////////////////////////////////////////////////////////

	DerivativeFunction::DerivativeFunction( const std::shared_ptr< IValue > & fixedTerm,
																					const std::shared_ptr< IValue > & derivativeTerm )
			: fixedTerm( fixedTerm ), derivativeTerm( derivativeTerm ) {}

//////////////////////////////////////////////////////////////////////////////
///  IElementLinear2D
//////////////////////////////////////////////////////////////////////////////

	IElementLinear2D::IElementLinear2D(
			const Node2D & t_Node1, const Node2D & t_Node2,
			const Node2D & t_Node3, const Node2D & t_Node4, const Material & t_Material ) :
			m_Material{ t_Material },
			m_Node{ { t_Node1, t_Node2, t_Node3, t_Node4 } },
			m_Global2D{ t_Node1, t_Node2, t_Node3, t_Node4 },
			m_QLECapacitance2D{ m_Global2D },
			m_QLEConductance2D{ m_Global2D } {
		auto matName = m_Material.name();
		auto & nodePool = NodePool::Instance();
		for( auto & node : m_Node ) {
			auto & poolNode = nodePool.getNode(node.getNodeNumber());
			poolNode.assignMaterial(matName);
		}
	}

	FenestrationCommon::SparceSquareMatrix< double > IElementLinear2D::conductanceMatrix() const {
		FenestrationCommon::SparceSquareMatrix< double > result{ numOfQuadrilateralNodes };
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

	FenestrationCommon::SparceSquareMatrix< double > IElementLinear2D::conductanceDerivativeMatrix() {
		FenestrationCommon::SparceSquareMatrix< double > result{ numOfQuadrilateralNodes };
		const auto numOfIntegrationPoints = IntegrationPoints2D::Instance().count2D();

		/// Integration matrix must be created every time
		std::vector< double > aDerivatives( numOfIntegrationPoints );

		auto count = 0u;
		m_QLEDerivativeConductance.clear();
		for ( const auto & cond : m_DerivativeConductance ) {
			for ( auto i = 0u; i < numOfIntegrationPoints; ++i ) {
				aDerivatives[ i ] = cond.derivativeTerm->value( m_Node[ i ].getState() );
			}
			m_QLEDerivativeConductance.emplace_back( m_Global2D );
			m_QLEDerivativeConductance[ count ].updateIntegrationMatrix( aDerivatives );
			++count;
		}

		/// Now rest of integration is performed as usual
		count = 0u;
		for ( const auto & cond : m_DerivativeConductance ) {
			std::vector< double > values( numOfIntegrationPoints );
			for ( auto i = 0u; i < numOfIntegrationPoints; ++i ) {
				values[ i ] = cond.fixedTerm->value( m_Node[ i ].getState() );
			}
			result += m_QLEDerivativeConductance[ count ].integrate( values );
			++count;
		}

		return result;
	}

	FenestrationCommon::SparceSquareMatrix< double > IElementLinear2D::capacitanceMatrix() const {
		FenestrationCommon::SparceSquareMatrix< double > result{ numOfQuadrilateralNodes };
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

	const Material & IElementLinear2D::getMaterial() const {
		return m_Material;
	}

	bool IElementLinear2D::haveBothNodes( const Node2D & t_Node1, const Node2D & t_Node2 ) const {
		bool node1Found = false;
		bool node2Found = false;
		for( auto & node : m_Node ) {
			node1Found = node1Found || node == t_Node1;
			node2Found = node2Found || node == t_Node2;
		}
		return node1Found && node2Found;
	}

	//////////////////////////////////////////////////////////////////////////////
	///  ElementThermalLinear2D
	//////////////////////////////////////////////////////////////////////////////

	ElementThermalLinear2D::ElementThermalLinear2D( const Node2D & t_Node1, const Node2D & t_Node2,
																									const Node2D & t_Node3, const Node2D & t_Node4,
																									const Material & mat ) :
			IElementLinear2D( t_Node1, t_Node2, t_Node3, t_Node4, mat ) {

		iValue waterFill = MaterialProperties::getWaterFill( mat );
		iValue airFill = MaterialProperties::getAirFill( mat );

		auto waterCapacitance = waterFill * ( Constants::Density_Water * Constants::Cp_Water );
		auto airCapacitance = airFill * ( Constants::Density_Air * Constants::Cp_Air );
		auto dryCapacitance = ( 1 - mat.porosity() ) * ( mat.density() * mat.heatCapacity() );

		auto capacitance = waterCapacitance + airCapacitance;
		capacitance = capacitance + dryCapacitance;

		m_Capacitance.push_back( capacitance );

		auto waterConductance = waterFill * Constants::K_Water;
		auto airConductance = airFill * Constants::K_Air;
		auto dryConductance = ( 1 - mat.porosity() ) * mat.thermalConductivity();

		auto conductance = waterConductance + airConductance;
		conductance = conductance + dryConductance;

		m_Conductance.push_back( conductance );

		/// Seems that it is not neccessary to create separate term when conductance coefficient
		/// is changing. Need to confirm that with rest of the team (Simon)
		///iValue one( std::make_shared< MoisThermFEM::Constant >( 1 ) );
		///m_DerivativeConductance.emplace_back( one, conductance );

	}

	//////////////////////////////////////////////////////////////////////////////
	///  ElementMoistureLinear2D
	//////////////////////////////////////////////////////////////////////////////

	ElementMoistureLinear2D::ElementMoistureLinear2D( const Node2D & t_Node1, const Node2D & t_Node2,
																										const Node2D & t_Node3, const Node2D & t_Node4,
																										const Material & mat ) :
			IElementLinear2D( t_Node1, t_Node2, t_Node3, t_Node4, mat ) {

		//////////////////////////////////////////////////////////////////////////////
		/// Creating conductance function for vapor
		//////////////////////////////////////////////////////////////////////////////
		iValue delta( std::make_shared< MoisThermFEM::Constant >( 2.5E-5 / mat.diffusionResistanceFactor() ) );
		iValue saturationFunction(
				std::make_shared< MoisThermFEM::SaturationFunction >( Property::temperature ) );

		m_Conductance.push_back( delta * saturationFunction );

		m_DerivativeConductance.emplace_back( delta, saturationFunction );

		//////////////////////////////////////////////////////////////////////////////
		/// Creating conductance function for liquid
		//////////////////////////////////////////////////////////////////////////////
		iValue suctionCurve(
				std::make_shared< MoisThermFEM::SuctionFunction >( mat.liquidTransportationCurve(),
																													 Property::humidity ) );
		m_Conductance.push_back( suctionCurve );

		//////////////////////////////////////////////////////////////////////////////
		/// Creating capacitance function
		//////////////////////////////////////////////////////////////////////////////
		iValue sorptionDerivative(
				std::make_shared< MoisThermFEM::TabularDerivative >( mat.sorptionCurve(),
																													 Property::humidity ) );

		m_Capacitance.push_back( sorptionDerivative );
	}
}
