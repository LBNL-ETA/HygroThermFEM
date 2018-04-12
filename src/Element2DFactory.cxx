#include "Element2DFactory.hxx"
#include "FEMunique.hxx"

namespace MoisThermFEM {

	void Element2DFactory::createThermalElement( const Node2D & t_Node1, const Node2D & t_Node2,
																							 const Node2D & t_Node3, const Node2D & t_Node4,
																							 const Material & mat ) {
		m_Elements.push_back(
				fem::make_unique< ElementThermalLinear2D >( t_Node1, t_Node2, t_Node3, t_Node4, mat ) );
	}

	void Element2DFactory::createMoistureElement( const Node2D & t_Node1, const Node2D & t_Node2,
																								const Node2D & t_Node3, const Node2D & t_Node4,
																								const Material & mat ) {
		m_Elements.push_back( fem::make_unique< ElementMoistureLinear2D >( t_Node1, t_Node2, t_Node3, t_Node4, mat ) );
		m_Linear = false;
	}
}
