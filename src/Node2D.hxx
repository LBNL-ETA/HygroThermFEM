#pragma once

#include <vector>
#include <set>
#include <string>
#include <initializer_list>

#include "State.hxx"
#include "Material.hxx"
#include "Node2D.hxx"
#include "Functions.hxx"


namespace HygroThermFEM
{
    ////////////////////////////////////////////////////////////////////////////
    ///   LocalPoint1D
    ////////////////////////////////////////////////////////////////////////////

    //! Holds local point coordinate in 1D
    struct LocalPoint1D
    {
        explicit LocalPoint1D(double t_ksi);

        LocalPoint1D(const LocalPoint1D & t_LocalPoint);

        double ksi{0};
    };

    ////////////////////////////////////////////////////////////////////////////
    ///   LocalPoint2D
    ////////////////////////////////////////////////////////////////////////////

    //! \brief Holds local point coordinates in 2D
    struct LocalPoint2D
    {
        LocalPoint2D(double t_ksi, double t_eta);

        LocalPoint2D(const LocalPoint2D & t_LocalPoint);

        double ksi{0};
        double eta{0};
    };

    //! \brief Main enumeration for accessing calculated or master properties.
    enum class Variable
    {
        temperature,
        //!< Temperature
        humidity,
        //!< Humidity
        pressure,
        //!< Pressure
        liquidPercent,
        //!< Percentage of water in liquid state.
        water,
        //!< Water content
        liquid,
        //!< Water content in liquid state.
        vapor,
        //!< Water content in vapor state.
        ice //!< Water content in frozen state.
    };

    //! \brief Enumerator for timestep count
    enum class Timestep
    {
        Current,
        //!< Current timestep
        Previous //!< Previous timestep
    };

    //! Simple Node2D interface
    class INode2D
    {
    public:
        INode2D() = delete;
        INode2D(std::size_t t_NodeNumber, double t_x, double t_y);
        virtual ~INode2D() = default;
        INode2D(const INode2D &) = default;
        [[nodiscard]] size_t getNodeNumber() const;
        double X() const;
        double Y() const;

        virtual void assignMaterial(const std::string & t_Material,
                                    double weightingCoefficient) = 0;

        virtual double property(Variable property,
                                Timestep iteration = Timestep::Current) const = 0;

        virtual void setStateProperty(BaseVariable t_Property,
                                      double t_value,
                                      bool updatePreviousValue = true) = 0;

        friend bool operator==(const INode2D & first, const INode2D & second);
        friend bool operator!=(const INode2D & first, const INode2D & second);

    protected:
        std::size_t m_NodeNumber{0};
        double m_x{0};
        double m_y{0};
    };

    ////////////////////////////////////////////////////////////////////////////
    ////   Node2D
    ////////////////////////////////////////////////////////////////////////////

    //! \brief Defines node in two dimensional coordinate space.
    //!
    //! Node is base construct of finite element model. Beside geometry description, node also
    //! contains information about state at its position and materials. Note that one node can have
    //! multiple materials assigned to it since one node can belong to multiple elements at same
    //! time.
    class Node2D : public INode2D
    {
    public:
        //! Base constructor of the node.
        Node2D(
            std::size_t t_NodeNumber,
            //!< Node index. Must be unique in entire finite element model
            double t_x,
            //!< Node's x coordinate
            double t_y,
            //!< Node's y coordinate
            const State & t_State //!< State of variables in the node
            );

        Node2D(const Node2D & t_Node) = default;
        Node2D(Node2D && t_Node) = default;
        ~Node2D() = default;
        Node2D & operator=(const Node2D & other) = default;
        Node2D & operator=(Node2D && other) = default;        

        //! Assigning material to current node.
        void assignMaterial(
            const std::string &
            t_Material,
            //!< SolidMaterial name. It must be assigned to MaterialPool first.
            double weightingCoefficient //!< Weighting coefficient that represents influence of the
                                        //!< material to current node.
            ) override;

        std::vector<std::string> getSolidMaterialNames() const;

        //! Returns back property of state variable. Variable can be from basic state (temperature,
        //! humidity or pressure) or water content
        double property(
            Variable property,
            //!< Variable for which value will be returned
            Timestep iteration = Timestep::Current //!< Timestep for which value will be returned
            ) const override;

        //! Sets the value of basic state property (temperature, humidity, pressure or liquid water
        //! percentage)
        void setStateProperty(
            BaseVariable t_Property,
            //!< Base state property for which value will be set
            double t_value,
            //!< New value that property will be set to.
            bool updatePreviousValue =
                true //!< Indicates wheter previous timestep value should be updated or not.
            ) override;

    private:
        //! Returns value of water content
        double waterContent(WaterContent content //!< Water content property
            ) const;

        //! Performs water content calculations and store it locally (This speeds up engine
        //! calculation)
        Water calcWaterContent();

        //! Update water content for every state (liquid, vapor and ice)
        void updateWaterContent();

        std::map<Timestep, State> m_State;

        struct MaterialContainer
        {
            MaterialContainer(const double weightingFactor, const IMaterial & material) :
                weightingFactor(weightingFactor),
                material(material)
            {}

            void addWeight(const double value)
            {
                weightingFactor += value;
            }

            double weightingFactor;
            const IMaterial & material;
        };

        /// Node can belong to multiple materials. This will be used to calculate secondary
        /// properties based on primary properties (water content depends on humidity)
        std::map<std::string, MaterialContainer> m_Materials;

        Water m_Water;
    };

    ////////////////////////////////////////////////////////////////////////////
    ////   INodes
    ////////////////////////////////////////////////////////////////////////////

    //! \brief Interface container that hold nodes in single place.
    //!
    //! In finite element model implemented in this engine, single element supports four nodes and
    //! single boundary condition support two nodes.
    //! This class is used to keep nodes in one place and also to provide some
    //! basic functionality to access necessary node data
    class INodes
    {
    public:
        INodes() = default;

        //! Construction of nodes storage from initializer list.
        INodes(std::initializer_list<std::reference_wrapper<INode2D>> t_Nodes);

        //! Returns properties for all nodes in the storage.
        [[nodiscard]] std::vector<double>
            properties(Variable property //!< Variable for which node values will be calculated.
                ) const;

        //! Simple operator[] overload for access to nodes storage by index.
        INode2D & operator[](std::size_t index //!< Node index
            ) const;

        //! Returns node indexes.
        [[nodiscard]] std::vector<std::size_t> getNodeIndexes() const;

        //! Returns number of nodes.
        [[nodiscard]] std::size_t size() const;

    protected:
        //! Construction of node storage for boundary conditions.
        INodes(Node2D & node1, Node2D & node2);

        //! Construction of node storage for element.
        INodes(Node2D & node1, Node2D & node2, Node2D & node3, Node2D & node4);

        std::vector<std::reference_wrapper<INode2D>> m_Nodes;
    };

    ////////////////////////////////////////////////////////////////////////////
    ////   LineNodes2D
    ////////////////////////////////////////////////////////////////////////////

    //! \brief Container for boundary condition.
    class LineNodes2D : public INodes
    {
    public:
        //! Container constructor
        LineNodes2D(Node2D & t_Node1, Node2D & t_Node2);
    };

    ////////////////////////////////////////////////////////////////////////////
    ////   QuadrilateralNodes2D
    ////////////////////////////////////////////////////////////////////////////

    //! \brief Container for quadrilateral element.
    class QuadrilateralNodes2D : public INodes
    {
    public:
        //! Container constructor
        QuadrilateralNodes2D(Node2D & t_Node1,
                             Node2D & t_Node2,
                             Node2D & t_Node3,
                             Node2D & t_Node4);
    };
} // namespace HygroThermFEM
