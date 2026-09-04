#include <functional>
#include <memory>
#include <vector>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "TestMaterials.hxx"

using HygroThermFEM::State;
using HygroThermFEM::Variable;

// Mass-conservation check for the moisture solver (defect D6).
//
// A closed strip (zero flux on every boundary -- no BCs applied) is initialised
// with a humidity gradient and left to relax. No moisture can enter or leave, so
// the total stored moisture integral(w dx dy) must stay constant. The old tangent
// capacity dw/dphi creates/destroys moisture where the sorption curve is steep;
// the mass-conservative secant capacity (SorptionSecantCapacity) keeps it constant.
namespace
{
    // 6 node columns x 1 element row, Cottaer Sandstone, moisture only.
    constexpr std::size_t nCols = 6;
    constexpr double dx = 0.02;
    constexpr double height = 0.05;

    std::size_t nodeIndex(std::size_t c, std::size_t r)
    {
        return c * 2 + r + 1;
    }

    // Total stored moisture: sum over elements of (mean nodal water) * element area.
    double totalMoisture(const std::vector<double> & water)
    {
        double total = 0.0;
        for(std::size_t c = 0; c < nCols - 1; ++c)
        {
            const double wMean = 0.25
                                 * (water[nodeIndex(c, 0) - 1] + water[nodeIndex(c + 1, 0) - 1]
                                    + water[nodeIndex(c + 1, 1) - 1] + water[nodeIndex(c, 1) - 1]);
            total += wMean * dx * height;
        }
        return total;
    }

    // Builds the closed strip with the given per-column humidity, relaxes it for
    // 24 hourly steps, and returns the largest relative drift of the moisture
    // integral.
    double closedStripMaxDrift(const std::function<double(std::size_t)> & initialHumidity)
    {
        HygroThermFEM::MultiDomain multiDomain(
          {.performThermal = false, .performMoisture = true});

        for(std::size_t c = 0; c < nCols; ++c)
        {
            const State state({.temperature = 20.0,
                               .humidity = initialHumidity(c),
                               .pressure = 101325.0,
                               .liquidPercent = 1.0});
            const double xPos = dx * static_cast<double>(c);
            multiDomain.nodes().createNode(
              {.index = nodeIndex(c, 0), .x = xPos, .y = 0.0, .state = state});
            multiDomain.nodes().createNode(
              {.index = nodeIndex(c, 1), .x = xPos, .y = height, .state = state});
        }

        const auto & material =
          multiDomain.materials().createSolidMaterial(TestHelper::CottaerSandstone());
        for(std::size_t c = 0; c < nCols - 1; ++c)
        {
            multiDomain.createElement({.node1 = nodeIndex(c, 0),
                                       .node2 = nodeIndex(c + 1, 0),
                                       .node3 = nodeIndex(c + 1, 1),
                                       .node4 = nodeIndex(c, 1),
                                       .material = material.name()});
        }

        // No boundary conditions -> zero flux on all edges (closed system).
        constexpr double dTime = 3600.0;
        constexpr unsigned nSteps = 24;

        auto humidities = multiDomain.nodes().properties(Variable::humidity);
        const double m0 = totalMoisture(multiDomain.nodes().properties(Variable::water));

        double maxDrift = 0.0;
        for(unsigned i = 0; i < nSteps; ++i)
        {
            humidities = multiDomain.moisture().transient(humidities, dTime).value().solution;
            // Advance the previous-timestep humidity as MultiDomain::transient does per step;
            // the secant capacity's conservation telescoping needs the true previous step.
            multiDomain.nodes().updateNodeHumidities(humidities, true);
            const double m = totalMoisture(multiDomain.nodes().properties(Variable::water));
            maxDrift = std::max(maxDrift, std::abs(m - m0) / m0);
        }

        return maxDrift;
    }
}   // namespace

TEST(MoistureMassConservation, ClosedStripConservesMoisture)
{
    // Humidity gradient 0.9 -> 0.6 across the columns; uniform temperature.
    const double maxDrift = closedStripMaxDrift(
      [](std::size_t c)
      { return 0.9 - 0.3 * static_cast<double>(c) / static_cast<double>(nCols - 1); });

    // With the mass-conservative secant capacity (exact nodal lumping) the closed-system
    // drift is machine precision (measured ~3e-14).
    EXPECT_LT(maxDrift, 1e-12);
}

TEST(MoistureMassConservation, ClosedStripConservesMoistureAboveLiquidData)
{
    // Humidity gradient 0.99 -> 0.96: water content 34.5-63 kg/m3, entirely above
    // Cottaer's first measured liquid transport point (w = 27). The drier gradient
    // above runs below that point, on a curve segment set by the zero-flooring
    // convention of the liquid table, so it cannot see the liquid coefficient at
    // all; this sibling exercises conservation where liquid transport is live.
    const double maxDrift = closedStripMaxDrift(
      [](std::size_t c)
      { return 0.99 - 0.03 * static_cast<double>(c) / static_cast<double>(nCols - 1); });

    EXPECT_LT(maxDrift, 1e-12);
}
