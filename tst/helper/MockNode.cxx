#include "MockNode.hxx"

// Mock classes for Node and Nodes so that we can test some functionality in real nodes.

namespace HygroThermFEM
{
    ////////////////////////////////////////////////////////////
    // StateValues
    ////////////////////////////////////////////////////////////
    StateValues::StateValues(const HygroThermFEM::Variable variable,
                             const double value,
                             const double prevTimestep) :
        variable(variable), value(value), prevTimestep(prevTimestep)
    {}

    ////////////////////////////////////////////////////////////
    // MockNode2D
    ////////////////////////////////////////////////////////////
    MockNode2D::MockNode2D(const size_t nodeNum, const double t_x, const double t_y) :
      INode2D(nodeNum, t_x, t_y),
      m_Property{{Variable::temperature, {{Timestep::Current, 0}, {Timestep::Previous, 0}}},
                 {Variable::humidity, {{Timestep::Current, 0}, {Timestep::Previous, 0}}}}
    {}

    MockNode2D::MockNode2D(const size_t nodeNum,
                           const double t_x,
                           const double t_y,
                           const std::vector<StateValues> & values) :
      INode2D(nodeNum, t_x, t_y)
    {
        for(const auto & val : values)
        {
            m_Property[val.variable] = {{Timestep::Current, val.value},
                                        {Timestep::Previous, val.prevTimestep}};
        }
    }

    MockNode2D::MockNode2D(const size_t nodeNum,
                           const double t_x,
                           const double t_y,
                           const StateValues & values) :
      INode2D(nodeNum, t_x, t_y),
      m_Property{
        {values.variable,
          {{Timestep::Current, values.value}, {Timestep::Previous, values.prevTimestep}}}}
    {}

    double MockNode2D::property(Variable variable, const Timestep timestep) const
    {
        return m_Property.at(variable).at(timestep);
    }

    void MockNode2D::assignMaterial(const IMaterial &, double)
    {}

    void MockNode2D::setStateProperty(BaseVariable, double, bool)
    {}

    ////////////////////////////////////////////////////////////////////////////
    ////   MockNodes2D
    ////////////////////////////////////////////////////////////////////////////
    MockNodes2D::MockNodes2D(MockNode2D & t_Node1, MockNode2D & t_Node2) :
      INodes{t_Node1, t_Node2}
    {}
}   // namespace HygroThermFEM