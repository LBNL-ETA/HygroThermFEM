#pragma once

#include <vector>
#include <set>
#include <string>
#include <initializer_list>

#include "State.hxx"
#include "Material.hxx"
#include "Node2D.hxx"

namespace MoisThermFEM
{
    ////////////////////////////////////////////////////////////////////////////
    ///   LocalPoint1D
    ////////////////////////////////////////////////////////////////////////////

    // Structure that holds data for one dimensional point in local coordinate system
    struct LocalPoint1D
    {
        explicit LocalPoint1D(const double t_ksi);

        LocalPoint1D(const LocalPoint1D & t_LocalPoint);

        double ksi{0};
    };

    ////////////////////////////////////////////////////////////////////////////
    ///   LocalPoint2D
    ////////////////////////////////////////////////////////////////////////////

    // Structure that holds data for two dimensional point in local coordinate system
    struct LocalPoint2D
    {
        LocalPoint2D(const double t_ksi, const double t_eta);

        LocalPoint2D(const LocalPoint2D & t_LocalPoint);

        double ksi{0};
        double eta{0};
    };

    enum class Property
	{
    	temperature,
    	humidity,
    	pressure,
    	liquidPercent,
    	water,
    	liquid,
    	vapor,
    	ice
	};

	enum class Timestep
	{
		Current,
		Previous
	};

    ////////////////////////////////////////////////////////////////////////////
    ////   Node2D
    ////////////////////////////////////////////////////////////////////////////

    // Defines nodal point in two dimensional cartesian space.
    class Node2D
    {
    public:
        Node2D(const std::size_t t_NodeNumber,
               const double t_x,
               const double t_y,
               const State & t_State);

        Node2D(const Node2D & t_Node) = default;
        Node2D & operator=(const Node2D & other) = default;
        friend bool operator==(const Node2D & first, const Node2D & second);
        friend bool operator!=(const Node2D & first, const Node2D & second);

        size_t getNodeNumber() const;

        double X() const;
        double Y() const;

        void assignMaterial( const std::string & t_Material, double weightingCoefficient );

        double property(const Property property, const Timestep iteration = Timestep::Current) const;
        double deltaProperty(const Property property) const;

		void setStateProperty( const StateProperty t_Property, double t_value );

    private:
		double waterContent( const WaterContent content ) const;

		double getStateProperty( const StateProperty t_Property,
								 const Timestep t_Iteration = Timestep::Current ) const;


        std::size_t m_NodeNumber{0};
        double m_x{0};
        double m_y{0};

        std::map<Timestep, State> m_State;

        /// Node can belong to multiple materials. This will be used to calculate secondary
        /// properties based on primary properties (water content depends on humidity)
        std::set<std::pair<double, std::reference_wrapper<const Material>>> m_Materials;
    };

    ////////////////////////////////////////////////////////////////////////////
    ////   INodes
    ////////////////////////////////////////////////////////////////////////////

    // Interface that holds all node data in single storage
    class INodes
    {
    public:
        INodes() = default;

        INodes(std::initializer_list<Node2D> t_Nodes);

        Node2D & getNode(const std::size_t Index);

        Node2D operator[](const std::size_t index) const;
        Node2D & operator[](const std::size_t index);

        std::vector<std::size_t> getNodeIndexes() const;

        std::size_t size() const;

    protected:
        std::vector<Node2D> m_Nodes;
    };

    ////////////////////////////////////////////////////////////////////////////
    ////   LineNodes2D
    ////////////////////////////////////////////////////////////////////////////

    // Class that store nodes which are part of some boundary conditions in 2D
    class LineNodes2D : public INodes
    {
    public:
        LineNodes2D(const Node2D & t_Node1, const Node2D & t_Node2);
    };

    ////////////////////////////////////////////////////////////////////////////
    ////   QuadrilateralNodes2D
    ////////////////////////////////////////////////////////////////////////////

    // Class that store nodal data which are part of elementsCreator
    class QuadrilateralNodes2D : public INodes
    {
    public:
        QuadrilateralNodes2D(const Node2D & t_Node1,
                             const Node2D & t_Node2,
                             const Node2D & t_Node3,
                             const Node2D & t_Node4);
    };

}   // namespace MoisThermFEM