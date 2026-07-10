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
}   // namespace

TEST(MoistureMassConservation, ClosedStripConservesMoisture)
{
    HygroThermFEM::MultiDomain multiDomain({.performThermal = false, .performMoisture = true});

    // Humidity gradient 0.9 -> 0.6 across the columns; uniform temperature.
    for(std::size_t c = 0; c < nCols; ++c)
    {
        const double phi = 0.9 - 0.3 * static_cast<double>(c) / static_cast<double>(nCols - 1);
        const State state({.temperature = 20.0,
                           .humidity = phi,
                           .pressure = 101325.0,
                           .liquidPercent = 1.0});
        multiDomain.nodes().createNode(
          {.index = nodeIndex(c, 0), .x = dx * static_cast<double>(c), .y = 0.0, .state = state});
        multiDomain.nodes().createNode(
          {.index = nodeIndex(c, 1), .x = dx * static_cast<double>(c), .y = height, .state = state});
    }

    const auto & material = multiDomain.materials().createSolidMaterial(TestHelper::CottaerSandstone());
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

    std::cout << "[MassConservation] initial moisture = " << m0
              << " kg/m, max relative drift = " << maxDrift << "\n";

    // With the mass-conservative secant capacity (exact nodal lumping) the closed-system
    // drift is machine precision (measured ~3e-14).
    EXPECT_LT(maxDrift, 1e-12);
}
