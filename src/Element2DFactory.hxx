#pragma once

#include "Elements2D.hxx"

namespace MoisThermFEM {

	class Node2D;
	class Material;

	class Element2DFactory : public ElementsLinear2D {
	public:
		void createThermalElement( const Node2D & t_Node1, const Node2D & t_Node2, const Node2D & t_Node3,
															 const Node2D & t_Node4, const Material & mat );
	};

}

