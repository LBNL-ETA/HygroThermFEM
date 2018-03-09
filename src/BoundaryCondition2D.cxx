#include <algorithm>
#include <functional>
#include <cmath>

#include "BoundaryCondition2D.hxx"
#include "Common.hxx"

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
	/// Flux BC
	////////////////////////////////////////////////////////

	FluxBC::FluxBC( Node2D & t_Node1, Node2D & t_Node2, const double t_Flux ) :
	IBCLinear2D( t_Node1, t_Node2 ), m_Flux( t_Flux )  {

	}

	FenestrationCommon::Vector< double > FluxBC::R_Vector() const {
		return m_PsiVector * m_Flux;
	}

	FenestrationCommon::SquareMatrix< double > FluxBC::H_Matrix() const {
		// Flux boundary conditions do not have H matrix (It is zero)
		return FenestrationCommon::SquareMatrix< double >( 4 );
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
										( std::pow( T, 2 ) + std::pow( m_RadiationTemperature, 2 ) ) *
										Constants::STEFANBOLTZMANN * m_Emissivity;
		}
		return result;
	}

	FenestrationCommon::Vector< double > BlackBodyRadiationBC::R_Vector() const {
		return m_PsiVector * HRadiative() * m_RadiationTemperature;
	}

	FenestrationCommon::SquareMatrix< double > BlackBodyRadiationBC::H_Matrix() const {
		return m_PsiPsiMatrix.mmultRows( HRadiative() );
	}

	/////////////////////////////////////////////////////
	/// MoistureBC
	/////////////////////////////////////////////////////

	MoistureBC::MoistureBC( const Node2D & t_Node1, const Node2D & t_Node2, const Material & material,
													const double t_ConvectiveCoefficient, const double t_AirHumidity,
													const double t_AirTemperature ) :
			IBCLinear2D( t_Node1, t_Node2 ), m_ConvectiveCoefficient( t_ConvectiveCoefficient ),
			m_AirHumidity( t_AirHumidity ), m_AirTemperature( t_AirTemperature ), m_Material( material ) {

	}

	FenestrationCommon::Vector< double > MoistureBC::R_Vector() const {
		using pValue = std::shared_ptr< IValue >;

		pValue saturation( std::make_shared< SaturationFunction >( Property::temperature ) );
		auto humidityCalculator = saturation * m_AirHumidity;

		humidityCalculator = humidityCalculator * m_Material.porosity();
		const auto humidityByVolume = humidityCalculator->value(
				State( m_AirTemperature, m_AirHumidity, 101325 ) );
		const auto coeff =
				m_ConvectiveCoefficient * humidityByVolume / ( Constants::Density_AIR * Constants::Cp_Air );
		return m_PsiVector * coeff;
	}

	FenestrationCommon::SquareMatrix< double > MoistureBC::H_Matrix() const {
		using pValue = std::shared_ptr< IValue >;

		pValue saturationFunction( std::make_shared< SaturationFunction >( Property::temperature ) );

		const auto humidityCoeff = m_Material.porosity() * saturationFunction;

		const auto humidityByVolume1 = humidityCoeff->value( m_Nodes[ 0 ].getState() );
		const auto humidityByVolume2 = humidityCoeff->value( m_Nodes[ 1 ].getState() );

		FenestrationCommon::Vector< double > coeffs {
			humidityByVolume1 * m_ConvectiveCoefficient / ( Constants::Density_AIR * Constants::Cp_Air ),
			humidityByVolume2 * m_ConvectiveCoefficient / ( Constants::Density_AIR * Constants::Cp_Air ),
		};

		return m_PsiPsiMatrix.mmultRows( coeffs );
	}
}