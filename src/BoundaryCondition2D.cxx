#include "BoundaryCondition2D.hxx"

namespace MoisThermFEM {

	////////////////////////////////////////////////////////
	// ConvectionBC
	////////////////////////////////////////////////////////

	ConvectionBC::ConvectionBC( const Node2D & t_Node1, const Node2D & t_Node2,
	                            const double t_ConvectionCoefficient, const double t_AirTemperature )
			: ILineLinear2D( t_Node1, t_Node2 ), m_ConvectionCoefficient( t_ConvectionCoefficient ),
			  m_AirTemperature( t_AirTemperature ) {

		const std::size_t numOfNodes = 2;
		// Create matrix A and vector R
		for( std::size_t i = 0; i < numOfIntegrationPoints(); ++i ) {
			for( std::size_t j = 0; j < numOfNodes; ++j ) {
				for( std::size_t k = 0; k < numOfNodes; ++k ) {
					m_matrixA[ j ][ k ] +=
							m_ConvectionCoefficient * m_Determinant * psi( i, j ) * psi( i, k );
				}
				m_Rvector[ j ] += m_ConvectionCoefficient * m_Determinant * m_AirTemperature * psi( i, j );
			}
		}

	}

	////////////////////////////////////////////////////////
	// TemperatureBC
	////////////////////////////////////////////////////////

	TemperatureBC::TemperatureBC( Node2D & t_Node1, Node2D & t_Node2,
	                              const double t_NodeTemperatures ) : ConvectionBC( t_Node1, t_Node2,
			1e18, t_NodeTemperatures ) {
		t_Node1.setProperty( Prop::temperature, t_NodeTemperatures );
		t_Node2.setProperty( Prop::temperature, t_NodeTemperatures );
	}

}