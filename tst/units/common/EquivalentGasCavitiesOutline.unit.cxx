#include <gtest/gtest.h>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>

#include "HygroThermFEM2D.hxx"
#include "GasCavities.hxx"

// The engine walks the outer edges of each frame cavity into one closed loop. A mesher can
// hand it a cavity that is no such thing: gmsh's quad recombination left one in two pieces
// joined at a node, and another with overlapping elements. The walk must then report the
// cavity. It used to read past the end of its edges and append until memory ran out.

namespace
{
    constexpr std::string_view cavityName{"Frame Cavity #7"};

    struct Model
    {
        HygroThermFEM::Nodes nodes;
        HygroThermFEM::Materials materials;
        HygroThermFEM::ElementsLinear2D elements;
    };

    void addUnitSquare(Model & model, const std::size_t firstNode, const double left)
    {
        model.nodes.createNode({.index = firstNode, .x = left, .y = 0.0, .state = {.temperature = 10.0}});
        model.nodes.createNode({.index = firstNode + 1u, .x = left + 1.0, .y = 0.0, .state = {.temperature = 10.0}});
        model.nodes.createNode({.index = firstNode + 2u, .x = left + 1.0, .y = 1.0, .state = {.temperature = 0.0}});
        model.nodes.createNode({.index = firstNode + 3u, .x = left, .y = 1.0, .state = {.temperature = 0.0}});
        model.elements.assignElement(std::make_unique<HygroThermFEM::ElementThermalLinear2D>(
          model.nodes,
          model.materials,
          firstNode,
          firstNode + 1u,
          firstNode + 2u,
          firstNode + 3u,
          std::string{cavityName},
          HygroThermFEM::PhysicsOptions{}));
    }

    void prepare(Model & model)
    {
        model.materials.createGas(std::string{cavityName}, HygroThermFEM::CavityStandard::ISO15099);
    }
}   // namespace

TEST(EquivalentGasCavitiesOutline, OneClosedOutlineBuildsTheCavity)
{
    Model model;
    prepare(model);
    addUnitSquare(model, 1u, 0.0);

    const HygroThermFEM::EquivalentGasCavities cavities(model.nodes, model.materials, model.elements);

    EXPECT_EQ(cavities.states().size(), 1u);
}

TEST(EquivalentGasCavitiesOutline, ACavityInTwoPiecesIsReportedByName)
{
    Model model;
    prepare(model);
    addUnitSquare(model, 1u, 0.0);
    addUnitSquare(model, 5u, 2.0);

    try
    {
        const HygroThermFEM::EquivalentGasCavities cavities(model.nodes, model.materials, model.elements);
        FAIL() << "two separate outlines must not build a cavity";
    }
    catch(const std::runtime_error & failure)
    {
        EXPECT_NE(std::string{failure.what()}.find(cavityName), std::string::npos) << failure.what();
    }
}

TEST(EquivalentGasCavitiesOutline, ACavityWithoutElementsIsReportedByName)
{
    Model model;
    prepare(model);

    EXPECT_THROW(HygroThermFEM::EquivalentGasCavities(model.nodes, model.materials, model.elements),
                 std::runtime_error);
}
