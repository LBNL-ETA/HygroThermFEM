#include <algorithm>
#include <functional>

#include "BoundaryCondition2D.hxx"

namespace MoisThermFEM {

	////////////////////////////////////////////////////////
	// ConvectionBC
	////////////////////////////////////////////////////////

	ConvectionBC::ConvectionBC( const Node2D & t_Node1, const Node2D & t_Node2,
															const double t_ConvectionCoefficient, const double t_AirTemperature )
			: IBCLinear2D( t_Node1, t_Node2 ), m_ConvectionCoefficient( t_ConvectionCoefficient ),
				m_AirTemperature( t_AirTemperature ) {

	}

	FenestrationCommon::Vector< double > ConvectionBC::R_Vector() const {
		return m_PsiVector * m_ConvectionCoefficient * m_AirTemperature;
	}

	FenestrationCommon::SquareMatrix< double > ConvectionBC::H_Matrix() const {
		return m_PsiPsiMatrix * m_ConvectionCoefficient;
	}

	////////////////////////////////////////////////////////
	/// TemperatureBC
	////////////////////////////////////////////////////////

	TemperatureBC::TemperatureBC( Node2D & t_Node1, Node2D & t_Node2,
																const double t_NodeTemperatures ) : ConvectionBC( t_Node1, t_Node2,
																																									1e18,
																																									t_NodeTemperatures ) {
		t_Node1.setProperty( Property::temperature, t_NodeTemperatures );
		t_Node2.setProperty( Property::temperature, t_NodeTemperatures );
	}

	TemperatureBC::TemperatureBC( Node2D & t_Node1, Node2D & t_Node2, const double t_Temp1,
																const double t_Temp2 ) : ConvectionBC( t_Node1, t_Node2, 1e18,
																																			 ( t_Temp1 + t_Temp2 ) / 2 ) {
		t_Node1.setProperty( Property::temperature, t_Temp1 );
		t_Node2.setProperty( Property::temperature, t_Temp2 );
	}

	////////////////////////////////////////////////////////
	/// BlackBodyRadiationBC
	////////////////////////////////////////////////////////

	BlackBodyRadiationBC::BlackBodyRadiationBC( const Node2D & t_Node1, const Node2D & t_Node2,
																							const double t_Emissivity,
																							const double t_RadiationTemperature )
			: IBCLinear2D( t_Node1, t_Node2, false ), m_RadiationTemperature { t_RadiationTemperature },
				m_Emissivity { t_Emissivity } {}

	FenestrationCommon::Vector< double > BlackBodyRadiationBC::HRadiative() const {
		FenestrationCommon::Vector< double > result( numOfBCNodes, 0 );
		for ( std::size_t j = 0; j < numOfBCNodes; ++j ) {
			double T = m_Nodes[ j ].getProperty( Property::temperature );
			result[ j ] = ( T + m_RadiationTemperature ) *
										( std::pow( T, 2 ) + std::pow( m_RadiationTemperature, 2 ) );
		}
		return result;
	}

	FenestrationCommon::Vector< double > BlackBodyRadiationBC::DHRadiative() const {
		FenestrationCommon::Vector< double > result( numOfBCNodes, 0 );
		for ( std::size_t j = 0; j < numOfBCNodes; ++j ) {
			double T = m_Nodes[ j ].getProperty( Property::temperature );
			result[ j ] = ( 3 * std::pow( T, 2 ) + 2 * m_RadiationTemperature * T +
											std::pow( m_RadiationTemperature, 2 ) );
		}
		return result;
	}

	FenestrationCommon::SquareMatrix< double > BlackBodyRadiationBC::D_HMatrix() {
		double integratedTemperature = getIntegratedProperty( Property::temperature );
		double integratedDeltaTemperature = getIntegratedDeltaProperty( Property::temperature );

		auto DhT = DHRadiative() * integratedTemperature;
		auto DhDt = DHRadiative() * integratedDeltaTemperature;
		auto VRadiation = DHRadiative() * m_RadiationTemperature;

		FenestrationCommon::SquareMatrix< double > C{ numOfBCNodes };
		FenestrationCommon::SquareMatrix< double > D{ numOfBCNodes };
		FenestrationCommon::SquareMatrix< double > E{ numOfBCNodes };

		for( auto i = 0u; i < numOfBCNodes; ++i ) {
			for( auto j = 0u; j < numOfBCNodes; ++j ) {
				C[ i ][ j ] = m_PsiPsiMatrix[ i ][ j ] * DhT[ i ];
				D[ i ][ j ] = m_PsiPsiMatrix[ i ][ j ] * DhDt[ i ];
				E[ i ][ j ] = m_PsiPsiMatrix[ i ][ j ] * VRadiation[ i ];
			}
		}

		return C + D + E;
	}
}