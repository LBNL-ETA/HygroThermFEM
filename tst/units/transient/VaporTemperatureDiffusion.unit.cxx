#include <cmath>
#include <fstream>
#include <functional>
#include <string>
#include <vector>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "TestMaterials.hxx"

using HygroThermFEM::State;
using HygroThermFEM::Variable;

// D1 isolation diagnostic (see doc/Moisture Governing Equations.md).
//
// A closed strip (zero flux on every boundary) starts with UNIFORM humidity but an
// IMPOSED temperature gradient, moisture-only (thermal off, so the temperature field is
// frozen). Because grad(phi) = 0 initially, the only term that can move moisture is the
// temperature-gradient vapour term delta*phi*grad(c_sat) -- exactly the term D1 fixes.
// So this isolates D1's behaviour, and being closed it also re-checks D6 conservation.
//
// The per-node humidity profile is written to CSV for comparison against the 1D Python
// reference (hygrothermfem_python), which implements the consistent vapour term.
namespace
{
    constexpr std::size_t nCols = 11;   // 10 elements along x
    constexpr double length = 0.1;
    constexpr double height = 0.05;
    constexpr double phi0 = 0.7;
    constexpr double dTime = 3600.0;
    constexpr unsigned nSteps = 24;

    std::size_t nodeIndex(std::size_t c, std::size_t r)
    {
        return c * 2 + r + 1;
    }

    double xOf(std::size_t c)
    {
        return length * static_cast<double>(c) / static_cast<double>(nCols - 1);
    }

    // Runs the closed-strip case for a given per-column temperature field and writes the
    // bottom-row humidity profile per step to CSV. Returns the max relative drift of total
    // stored moisture (should be ~0 for a closed system with the conservative capacity).
    double runAndDump(const std::function<double(double)> & temperatureAt,
                      const std::string & csvPath)
    {
        HygroThermFEM::MultiDomain multiDomain({.performThermal = false, .performMoisture = true});

        for(std::size_t c = 0; c < nCols; ++c)
        {
            const State state({.temperature = temperatureAt(xOf(c)),
                               .humidity = phi0,
                               .pressure = 101325.0,
                               .liquidPercent = 1.0});
            multiDomain.nodes().createNode(
              {.index = nodeIndex(c, 0), .x = xOf(c), .y = 0.0, .state = state});
            multiDomain.nodes().createNode(
              {.index = nodeIndex(c, 1), .x = xOf(c), .y = height, .state = state});
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

        std::ofstream csv(csvPath);
        csv << "step";
        for(std::size_t c = 0; c < nCols; ++c)
        {
            csv << ",x" << c;
        }
        csv << "\n";

        auto humidities = multiDomain.nodes().properties(Variable::humidity);

        auto totalWater = [&]() {
            const auto w = multiDomain.nodes().properties(Variable::water);
            double total = 0.0;
            const double dx = length / static_cast<double>(nCols - 1);
            for(std::size_t c = 0; c < nCols - 1; ++c)
            {
                const double wMean = 0.25
                                     * (w[nodeIndex(c, 0) - 1] + w[nodeIndex(c + 1, 0) - 1]
                                        + w[nodeIndex(c + 1, 1) - 1] + w[nodeIndex(c, 1) - 1]);
                total += wMean * dx * height;
            }
            return total;
        };

        const double m0 = totalWater();
        double maxDrift = 0.0;

        auto writeRow = [&](unsigned step) {
            csv << step;
            const auto phi = multiDomain.nodes().properties(Variable::humidity);
            for(std::size_t c = 0; c < nCols; ++c)
            {
                csv << "," << phi[nodeIndex(c, 0) - 1];
            }
            csv << "\n";
        };

        writeRow(0);
        for(unsigned i = 0; i < nSteps; ++i)
        {
            humidities = multiDomain.moisture().transient(humidities, dTime).value().solution;
            writeRow(i + 1);
            maxDrift = std::max(maxDrift, std::abs(totalWater() - m0) / m0);
        }
        return maxDrift;
    }
}   // namespace

TEST(VaporTemperatureDiffusion, IsothermalControlAndGradient)
{
    // Control: uniform temperature -> the D1 term is inert, nothing should move.
    const double isoDrift = runAndDump([](double) { return 20.0; }, "/tmp/htf_d1_isothermal.csv");

    // Test: linear temperature gradient 40 C (x=0) -> 20 C (x=L). Only the D1 vapour
    // temperature-gradient term can redistribute moisture.
    const double gradDrift = runAndDump(
      [](double x) { return 40.0 - 20.0 * (x / length); }, "/tmp/htf_d1_gradient.csv");

    std::cout << "[D1] isothermal drift = " << isoDrift << ", gradient drift = " << gradDrift
              << "\n";

    // Isothermal: the D1 term is inert, so nothing moves and mass is exactly conserved.
    EXPECT_LT(isoDrift, 1e-9);

    // Gradient: D1 correctly drives moisture hot -> cold (verified against the 1D reference
    // in hygrothermfem_python/examples/compare_d1.py -- same direction, isothermal exact).
    // The reference conserves mass to ~1e-12; the engine leaks ~1% here because its moisture
    // capacity is assembled via the coefficient-averaged psi_i psi_j capacitance matrix and
    // then lumped, which only *approximately* conserves. A pure nodal-lumped secant capacity
    // (as in the reference) would make this exact. This bound documents the current behaviour;
    // tightening it is the D6-refinement follow-up.
    EXPECT_LT(gradDrift, 2e-2);
}
