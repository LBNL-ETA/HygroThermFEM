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

	std::vector< double > ConvectionBC::R_Vector() const {
		std::vector< double > A;
		std::transform( m_PsiVector.begin(), m_PsiVector.end(), std::back_inserter( A ),
										std::bind1st( std::multiplies< double >(),
																	m_ConvectionCoefficient *
																	m_AirTemperature ) );
		return A;
	}

	FenestrationCommon::SquareMatrix< double > ConvectionBC::H_Matrix() const {
		return m_PsiPsiMatrix.mult( m_ConvectionCoefficient );
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
																							const double t_Emissivity )
			: IBCLinear2D( t_Node1,
										 t_Node2,
										 true ), m_Emissivity( t_Emissivity ) {}
}