#include <gtest/gtest.h>
#include <memory>
#include <optional>
#include <set>

#include "HygroThermFEM2D.hxx"
#include "TestMaterials.hxx"

// The frame cavity of ModelWithFrameCavity1, solved with and without emissivity overrides on
// its boundary segments. Without an override every segment takes the emissivity of the
// surrounding solid material, which is the solid material default of 0.9. Overriding every
// segment to 0.9 must therefore reproduce the no-override result exactly, and overriding to
// 0.0 removes the radiative part and must lower the effective conductivity.
//
// The steady-state cases add a second convection boundary on the right edge so the cavity
// sits in a temperature gradient, and check that the steady coupling loop updates the cavity
// at all (its conductivity leaves the gas default) and honours the overrides the same way.
// The reported cavity state must describe the same evaluation.

namespace
{
    using HygroThermFEM::State;

    constexpr std::size_t gridSize{5u};
    const std::set<std::size_t> cavityElements{6u, 7u, 10u};
    constexpr double solidDefaultEmissivity{0.9};
    const std::string cavityName{"Frame Cavity 1"};

    enum class Solve
    {
        Transient,
        Steady
    };

    struct ElementNodes
    {
        std::size_t node1;
        std::size_t node2;
        std::size_t node3;
        std::size_t node4;
    };

    struct RunResult
    {
        double conductivity{0.0};
        std::vector<HygroThermFEM::GasCavityState> states;
    };

    ElementNodes elementNodes(const std::size_t xIndex, const std::size_t yIndex)
    {
        return {.node1 = xIndex * gridSize + yIndex - gridSize,
                .node2 = xIndex * gridSize + yIndex - gridSize + 1u,
                .node3 = xIndex * gridSize + (yIndex + 1u),
                .node4 = xIndex * gridSize + yIndex};
    }

    void createGridNodes(HygroThermFEM::MultiDomain & multiDomain)
    {
        const std::vector<double> gridX{0, 0.01, 0.02, 0.03, 0.04};
        const std::vector<double> gridY{0, 0.05, 0.1, 0.15, 0.2};
        const State state({.temperature = 20.0, .humidity = 0.0, .pressure = 101325.0});
        for(const auto yVal : gridY)
        {
            for(const auto xVal : gridX)
            {
                multiDomain.nodes().createNode({.x = xVal, .y = yVal, .state = state});
            }
        }
    }

    void overrideElementEdges(HygroThermFEM::MultiDomain & multiDomain,
                              const ElementNodes & nodes,
                              const double emissivity)
    {
        multiDomain.setCavitySurfaceEmissivity(nodes.node1, nodes.node2, emissivity);
        multiDomain.setCavitySurfaceEmissivity(nodes.node2, nodes.node3, emissivity);
        multiDomain.setCavitySurfaceEmissivity(nodes.node3, nodes.node4, emissivity);
        multiDomain.setCavitySurfaceEmissivity(nodes.node4, nodes.node1, emissivity);
    }

    void createElements(HygroThermFEM::MultiDomain & multiDomain,
                        const std::string & solidName,
                        const std::optional<double> surfaceEmissivity)
    {
        std::size_t elementNumber{0u};
        for(std::size_t xIndex = 1u; xIndex < gridSize; ++xIndex)
        {
            for(std::size_t yIndex = 1u; yIndex < gridSize; ++yIndex)
            {
                ++elementNumber;
                const auto nodes{elementNodes(xIndex, yIndex)};
                const bool isCavity{cavityElements.contains(elementNumber)};
                multiDomain.createElement({.node1 = nodes.node1,
                                           .node2 = nodes.node2,
                                           .node3 = nodes.node3,
                                           .node4 = nodes.node4,
                                           .material = isCavity ? cavityName : solidName});
                if(isCavity && surfaceEmissivity.has_value())
                {
                    overrideElementEdges(multiDomain, nodes, surfaceEmissivity.value());
                }
            }
        }
    }

    //! Convection boundary along one vertical grid column, from the bottom node up.
    void createColumnConvection(HygroThermFEM::MultiDomain & multiDomain,
                                const std::size_t column,
                                const double airTemperature)
    {
        const HygroThermFEM::FixedBCHCCoefficients bcCoeff{airTemperature, 30.0};
        std::vector<std::size_t> bcNodes;
        for(std::size_t idx = 0u; idx < gridSize; ++idx)
        {
            bcNodes.push_back(idx * gridSize + column);
        }
        for(std::size_t idx = 1u; idx < bcNodes.size(); ++idx)
        {
            multiDomain.thermal().createBC_FixedHc(bcNodes[idx - 1u], bcNodes[idx], bcCoeff);
        }
    }

    //! Runs the cavity model of ModelWithFrameCavity1, with every cavity boundary segment
    //! overridden to surfaceEmissivity when given. The transient run is that test's; the steady
    //! run adds a warm right edge so the cavity sees a gradient.
    RunResult runCavityModel(const std::optional<double> surfaceEmissivity, const Solve solve)
    {
        HygroThermFEM::MultiDomain multiDomain({.performThermal = true, .performMoisture = false});
        createGridNodes(multiDomain);

        auto params = TestHelper::CottaerSandstone();
        params.name = "Material 1";
        params.liquidTransportCurve = {{0, 0}, {180, 2e-6}};
        params.sorptionCurve = {{0, 0}, {1, 180}};
        const auto & solidMaterial = multiDomain.materials().createSolidMaterial(params);
        auto & frameCavity =
          multiDomain.materials().createGas(cavityName, HygroThermFEM::CavityStandard::ISO15099);

        createElements(multiDomain, solidMaterial.name(), surfaceEmissivity);
        createColumnConvection(multiDomain, 1u, 0.0);

        if(solve == Solve::Transient)
        {
            constexpr auto dTime = 360;
            constexpr auto nSteps = 10;
            multiDomain.thermal().transientMultiStep(HygroThermFEM::Variable::temperature, dTime, nSteps);
        }
        else
        {
            createColumnConvection(multiDomain, gridSize, 20.0);
            multiDomain.steadyState();
        }

        return {.conductivity = frameCavity.thermalConductivityMoistureAndTemperatureDependent().getCurve()[0].y,
                .states = multiDomain.gasCavityStates()};
    }

    double cavityConductivity(const std::optional<double> surfaceEmissivity, const Solve solve)
    {
        return runCavityModel(surfaceEmissivity, solve).conductivity;
    }
}   // namespace

TEST(TestModelWithFrameCavitySurfaceEmissivity, NoOverrideMatchesModelWithFrameCavity1)
{
    SCOPED_TRACE("Begin Test: Frame cavity without surface emissivity overrides.");

    const auto conductivity{cavityConductivity(std::nullopt, Solve::Transient)};

    EXPECT_NEAR(0.098654, conductivity, 1e-6);
}

TEST(TestModelWithFrameCavitySurfaceEmissivity, MaterialEmissivityOverrideEqualsNoOverride)
{
    SCOPED_TRACE("Begin Test: Overriding every segment to the solid material's own emissivity.");

    const auto withoutOverride{cavityConductivity(std::nullopt, Solve::Transient)};
    const auto withSameOverride{cavityConductivity(solidDefaultEmissivity, Solve::Transient)};

    EXPECT_DOUBLE_EQ(withoutOverride, withSameOverride);
}

TEST(TestModelWithFrameCavitySurfaceEmissivity, ZeroEmissivityLowersConductivity)
{
    SCOPED_TRACE("Begin Test: Overriding every segment to 0.0 removes radiation from the cavity.");

    const auto withoutOverride{cavityConductivity(std::nullopt, Solve::Transient)};
    const auto withZeroOverride{cavityConductivity(0.0, Solve::Transient)};

    EXPECT_LT(withZeroOverride, withoutOverride);
}

TEST(TestModelWithFrameCavitySurfaceEmissivity, SteadyStateUpdatesCavity)
{
    SCOPED_TRACE("Begin Test: The steady coupling loop recomputes the cavity conductivity.");

    constexpr auto gasDefaultConductivity{0.05};
    const auto conductivity{cavityConductivity(std::nullopt, Solve::Steady)};

    EXPECT_NE(gasDefaultConductivity, conductivity);
}

TEST(TestModelWithFrameCavitySurfaceEmissivity, SteadyStateHonoursOverride)
{
    SCOPED_TRACE("Begin Test: Steady state with overrides 0.9 and 0.0 on every segment.");

    const auto withoutOverride{cavityConductivity(std::nullopt, Solve::Steady)};
    const auto withSameOverride{cavityConductivity(solidDefaultEmissivity, Solve::Steady)};
    const auto withZeroOverride{cavityConductivity(0.0, Solve::Steady)};

    EXPECT_DOUBLE_EQ(withoutOverride, withSameOverride);
    EXPECT_LT(withZeroOverride, withoutOverride);
}

TEST(TestModelWithFrameCavitySurfaceEmissivity, SteadyStateReportsCavityState)
{
    SCOPED_TRACE("Begin Test: The reported cavity state matches the evaluation of the cavity.");

    constexpr double overriddenEmissivity{0.345};
    const auto result{runCavityModel(overriddenEmissivity, Solve::Steady)};

    ASSERT_EQ(1u, result.states.size());
    const auto & state{result.states.front()};
    EXPECT_EQ(cavityName, state.gasName);
    EXPECT_DOUBLE_EQ(result.conductivity, state.conductivity);
    EXPECT_GT(state.nusselt, 0.0);
    EXPECT_NEAR(overriddenEmissivity, state.side1.emissivity, 1e-12);
    EXPECT_NEAR(overriddenEmissivity, state.side2.emissivity, 1e-12);
    // Cold air on the left (0 C) and warm air on the right (20 C): the warm side is the right
    // one, so the screen heat flow runs to the left, and the sides are reported in Kelvin.
    EXPECT_EQ(KeffCavity::ScreenFlow::Left, state.heatFlow);
    EXPECT_GT(state.side1.temperature, state.side2.temperature);
    EXPECT_GT(state.side2.temperature, 270.0);
}
