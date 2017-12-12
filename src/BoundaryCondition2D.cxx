#include "BoundaryCondition2D.hxx"
#include "NodePool.hxx"

namespace MoisThermFEM {
  
  ////////////////////////////////////////////////////////
  // ConvectionBC
  ////////////////////////////////////////////////////////

  ConvectionBC::ConvectionBC( 
    Node2D const & t_Node1, 
    Node2D const & t_Node2, 
    double const t_ConvectionCoefficient, 
    double const t_AirTemperature ) :
    ILineLinear2D( t_Node1, t_Node2 ), 
    m_ConvectionCoefficient( t_ConvectionCoefficient ),
    m_AirTemperature( t_AirTemperature ) {

    // Create matrix A and vector R
    for( std::size_t i = 0; i < numOfIntegrationPoints(); ++i ) {
      for( std::size_t j = 0; j < 2; ++j ) {
        for( std::size_t k = 0; k < 2; ++k ) {
          m_matrixA[ j ][ k ] += m_ConvectionCoefficient * m_Determinant * Psi( i, j ) * Psi( i, k );
        }
        m_Rvector[ j ] += m_ConvectionCoefficient * m_Determinant * m_AirTemperature * Psi( i, j );
      }
    }

  }

  ////////////////////////////////////////////////////////
  // TemperatureBC
  ////////////////////////////////////////////////////////

  TemperatureBC::TemperatureBC( 
    Node2D & t_Node1, 
    Node2D & t_Node2, 
    const double t_NodeTemperatures ) : ConvectionBC( t_Node1, t_Node2, 1e18, t_NodeTemperatures ) {
    t_Node1.temperature = t_NodeTemperatures;
    t_Node2.temperature = t_NodeTemperatures;
  }

}