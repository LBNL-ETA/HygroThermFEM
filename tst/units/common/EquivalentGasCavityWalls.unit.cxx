#include <gtest/gtest.h>
#include <cstddef>
#include <memory>
#include <string>

#include "HygroThermFEM2D.hxx"
#include "GasCavities.hxx"

// A frame cavity side takes its emissivity from the walls on it. A wall is an edge of a solid
// element. An outline segment with no solid element on it opens into a neighbouring cavity (or
// out of the model): radiation passes through, so it counts as black, emissivity 1. That is the
// value THERM gives Conrad for an opening between two frame cavities, and the reference here.
//
// The second test is THERM's sample-sill.thmz in miniature: two cavities side by side between
// two rows of wood. Both end nodes of the opening between them belong to wood elements, yet no
// wood lies on the opening itself. The engine used to look for a material the two end nodes
// share, found the wood, and gave the opening the wood's emissivity.

namespace
{
    constexpr double woodEmissivity{0.6};
    constexpr double openingEmissivity{1.0};
    const std::string wood{"Wood"};
    const std::string leftCavity{"Left cavity"};
    const std::string rightCavity{"Right cavity"};

    struct Model
    {
        HygroThermFEM::Nodes nodes;
        HygroThermFEM::Materials materials;
        HygroThermFEM::ElementsLinear2D elements;
    };

    //! Nodes on a grid with unit spacing, numbered row by row from 1. The temperature falls from
    //! left to right, so the heat in every cavity flows to the right.
    void addGrid(Model & model, const std::size_t columns, const std::size_t rows)
    {
        for(std::size_t row = 0u; row <= rows; ++row)
        {
            for(std::size_t column = 0u; column <= columns; ++column)
            {
                const auto xCoordinate{static_cast<double>(column)};
                model.nodes.createNode({.index = 1u + row * (columns + 1u) + column,
                                        .x = xCoordinate,
                                        .y = static_cast<double>(row),
                                        .state = {.temperature = 20.0 - 5.0 * xCoordinate}});
            }
        }
    }

    //! The unit square whose lower left corner is grid point (column, row).
    void addSquare(Model & model, const std::size_t columns, const std::size_t column, const std::size_t row, const std::string & material)
    {
        const std::size_t lowerLeft{1u + row * (columns + 1u) + column};
        const std::size_t upperLeft{lowerLeft + columns + 1u};
        model.elements.assignElement(std::make_unique<HygroThermFEM::ElementThermalLinear2D>(
          model.nodes, model.materials, lowerLeft, lowerLeft + 1u, upperLeft + 1u, upperLeft, material, HygroThermFEM::PhysicsOptions{}));
    }

    void createWood(Model & model)
    {
        // The element reads the sorption curve when it is built; a dry solid holds no water.
        model.materials.createSolidMaterial({.name = wood,
                                             .thermalConductivityDry = 0.14,
                                             .sorptionCurve = {{0.0, 0.0}, {1.0, 0.0}},
                                             .emissivity = woodEmissivity});
    }

    void createCavity(Model & model, const std::string & name)
    {
        model.materials.createGas(name, HygroThermFEM::CavityStandard::ISO15099);
    }
}   // namespace

TEST(EquivalentGasCavityWalls, AWallTakesTheEmissivityOfTheSolidOnIt)
{
    // | wood | cavity |   heat flows to the right; nothing lies right of the cavity.
    Model model;
    createWood(model);
    createCavity(model, leftCavity);
    addGrid(model, 2u, 1u);
    addSquare(model, 2u, 0u, 0u, wood);
    addSquare(model, 2u, 1u, 0u, leftCavity);

    HygroThermFEM::EquivalentGasCavities cavities(model.nodes, model.materials, model.elements);
    cavities.update();
    const auto states{cavities.states()};

    ASSERT_EQ(states.size(), 1u);
    EXPECT_EQ(states[0].heatFlow, KeffCavity::ScreenFlow::Right);
    EXPECT_NEAR(states[0].side1.emissivity, woodEmissivity, 1e-12);
    EXPECT_NEAR(states[0].side2.emissivity, openingEmissivity, 1e-12);
}

TEST(EquivalentGasCavityWalls, AnOpeningBetweenTwoCavitiesIsBlackEvenWithWoodAtBothEnds)
{
    // | wood   | wood    |
    // | left   | right   |   cavities; the opening between them is x = 1.
    // | wood   | wood    |
    Model model;
    createWood(model);
    createCavity(model, leftCavity);
    createCavity(model, rightCavity);
    addGrid(model, 2u, 3u);
    addSquare(model, 2u, 0u, 0u, wood);
    addSquare(model, 2u, 1u, 0u, wood);
    addSquare(model, 2u, 0u, 1u, leftCavity);
    addSquare(model, 2u, 1u, 1u, rightCavity);
    addSquare(model, 2u, 0u, 2u, wood);
    addSquare(model, 2u, 1u, 2u, wood);

    HygroThermFEM::EquivalentGasCavities cavities(model.nodes, model.materials, model.elements);
    cavities.update();
    const auto states{cavities.states()};

    ASSERT_EQ(states.size(), 2u);
    for(const auto & state : states)
    {
        EXPECT_EQ(state.heatFlow, KeffCavity::ScreenFlow::Right) << state.gasName;
        EXPECT_NEAR(state.side1.emissivity, openingEmissivity, 1e-12) << state.gasName;
        EXPECT_NEAR(state.side2.emissivity, openingEmissivity, 1e-12) << state.gasName;
    }
}
