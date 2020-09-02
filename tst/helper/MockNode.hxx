#pragma once

#include "HygroThermFEM2D.hxx"

// Mock object that will be used to test tables
namespace HygroThermFEM
{
    ////////////////////////////////////////////////////////////
    // StateValues
    ////////////////////////////////////////////////////////////
    struct StateValues
    {
        StateValues(Variable variable, double value, double prevTimestep = 0);

        Variable variable{Variable::temperature};
        double value{0};
        double prevTimestep{0};
    };

    ////////////////////////////////////////////////////////////
    // MockNode2D
    ////////////////////////////////////////////////////////////
    class MockNode2D final : public INode2D
    {
    public:
        explicit MockNode2D(size_t nodeNum = 0, double t_x = 0, double t_y = 0);

        MockNode2D(size_t nodeNum, double t_x, double t_y, const std::vector<StateValues> & values);

        ~MockNode2D() override = default;
        MockNode2D(const MockNode2D & mockNode) = default;
        MockNode2D(MockNode2D && mockNode) = default;
        MockNode2D & operator=(const MockNode2D &) = default;
        MockNode2D & operator=(MockNode2D &&) = default;

        void assignMaterial(const std::string &, double) override;

        void setStateProperty(BaseVariable, double, bool) override;

        explicit MockNode2D(size_t nodeNum, double t_x, double t_y, const StateValues & values);

        [[nodiscard]] double property(Variable variable, Timestep timestep) const override;

    private:
        std::map<Variable, std::map<Timestep, double>> m_Property;
    };

    ////////////////////////////////////////////////////////////////////////////
    ////   MockNodes2D
    ////////////////////////////////////////////////////////////////////////////

    class MockNodes2D : public INodes
    {
    public:
        //! Container constructor
        MockNodes2D(MockNode2D & t_Node1, MockNode2D & t_Node2);
    };
}   // namespace HygroThermFEM