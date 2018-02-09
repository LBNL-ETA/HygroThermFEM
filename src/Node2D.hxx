#pragma once

#include <vector>
#include <initializer_list>
#include <map>

#include "State.hxx"

namespace MoisThermFEM {

	////////////////////////////////////////////////////////////////////////////
	//   LocalPoint1D
	////////////////////////////////////////////////////////////////////////////

	// Structure that holds data for one dimensional point in local coordinate system
	struct LocalPoint1D {
		explicit LocalPoint1D( const double t_ksi );

		LocalPoint1D( const LocalPoint1D & t_LocalPoint );

		double ksi { 0 };

	};

	////////////////////////////////////////////////////////////////////////////
	//   LocalPoint2D
	////////////////////////////////////////////////////////////////////////////

	// Structure that holds data for two dimensional point in local coordinate system
	struct LocalPoint2D {
		LocalPoint2D( const double t_ksi, const double t_eta );

		LocalPoint2D( const LocalPoint2D & t_LocalPoint );

		double ksi { 0 };
		double eta { 0 };

	};

	////////////////////////////////////////////////////////////////////////////
	///   Node2D
	////////////////////////////////////////////////////////////////////////////

	// Defines nodal point in two dimensional cartesian space.
	class Node2D {
	public:
		Node2D( const std::size_t t_NodeNumber, const double t_x, const double t_y,
						const State & t_State );

		Node2D( const Node2D & t_Node );
		Node2D & operator=( const Node2D & other );

		size_t getNodeNumber() const;

		double X() const;
		double Y() const;

		double getProperty( const Property t_Property,
												const Iteration t_Iteration = Iteration::Current ) const;
		void setProperty( const Property t_Property, double t_value );
		double getDeltaProperty( const Property t_Property ) const;

	private:
		std::size_t m_NodeNumber { 0 };
		double m_x { 0 };
		double m_y { 0 };

		State m_State;

	};

	////////////////////////////////////////////////////////////////////////////
	//   INodesStorage
	////////////////////////////////////////////////////////////////////////////

	// Interface that holds all node data in single storage
	class INodesStorage {
	public:
		INodesStorage() = default;

		INodesStorage( std::initializer_list< Node2D > t_Nodes );

		Node2D getNode( const std::size_t Index ) const;

		Node2D operator[]( const std::size_t index ) const;

		std::vector< std::size_t > getNodeIndexes() const;

	protected:
		std::vector< Node2D > m_Nodes;

	};

	////////////////////////////////////////////////////////////////////////////
	//   LineNodes2D
	////////////////////////////////////////////////////////////////////////////

	// Class that store nodes which are part of some boundary conditions in 2D
	class LineNodes2D : public INodesStorage {
	public:
		LineNodes2D( const Node2D & t_Node1, const Node2D & t_Node2 );

	};

	////////////////////////////////////////////////////////////////////////////
	//   QuadrilateralNodes2D
	////////////////////////////////////////////////////////////////////////////

	// Class that store nodal data which are part of elements
	class QuadrilateralNodes2D : public INodesStorage {
	public:
		QuadrilateralNodes2D( const Node2D & t_Node1, const Node2D & t_Node2, const Node2D & t_Node3,
													const Node2D & t_Node4 );
	};

}