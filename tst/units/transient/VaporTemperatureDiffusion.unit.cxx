#include <array>
#include <cmath>
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
// The full-profile comparison against the 1D reference solver lives in the validation
// book (hygrothermfem_python, cases vapor_isothermal_cottaer and vapor_gradient_cottaer);
// what is asserted here is the closed system's own invariant, its stored moisture.
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

    struct StripRun
    {
        //! Max relative drift of total stored moisture over the run (should be ~0 for a
        //! closed system with the conservative capacity).
        double maxDrift;
        //! Final bottom-row humidity profile, one value per column.
        std::vector<double> finalBottomRow;
    };

    // Runs the closed-strip case for a given per-column temperature field.
    StripRun runClosedStrip(const std::function<double(double)> & temperatureAt)
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

        for(unsigned i = 0; i < nSteps; ++i)
        {
            humidities = multiDomain.moisture().transient(humidities, dTime).value().solution;
            // Advance the nodes' previous-timestep humidity, exactly as MultiDomain::transient
            // does after each completed step. Without this the secant storage capacity keeps
            // measuring against the INITIAL state instead of the previous step, which breaks
            // its mass-conservation telescoping for every step after the first.
            multiDomain.nodes().updateNodeHumidities(humidities, true);
            maxDrift = std::max(maxDrift, std::abs(totalWater() - m0) / m0);
        }

        std::vector<double> bottomRow;
        bottomRow.reserve(nCols);
        for(std::size_t col = 0; col < nCols; ++col)
        {
            bottomRow.push_back(humidities[nodeIndex(col, 0) - 1]);
        }
        return {.maxDrift = maxDrift, .finalBottomRow = std::move(bottomRow)};
    }
}   // namespace

TEST(VaporTemperatureDiffusion, GradientVaporOnly)
{
    // Same closed strip and imposed gradient, but with liquid transportation excluded.
    // With the temperature frozen the vapour coefficients are constant during the
    // nonlinear iteration, so any remaining moisture drift here is an operator-level
    // conservation error, not iteration (Picard) quality -- this splits the full-physics
    // gradient drift into its liquid/iteration and operator parts.
    HygroThermFEM::SimulationProperties::Instance().setCalculationParameters(
      true, false, false, false, false);

    const double vaporOnlyDrift =
      runClosedStrip([](double x) { return 40.0 - 20.0 * (x / length); }).maxDrift;

    HygroThermFEM::SimulationProperties::Instance().resetCalculationParameters();

    std::cout << "[D1] vapor-only gradient drift = " << vaporOnlyDrift << "\n";

    // Bounded by nonlinear-iteration residual only (measured 2.3e-8 over 24 steps).
    // NOTE: for Cottaer this case is numerically identical to the full-physics gradient
    // case, because the logarithmically interpolated liquid transport curve is exactly
    // zero across its first segment (anchor point (0, 0)) -- there is no liquid
    // transport at all below w = 27 kg/m^3 (phi ~ 0.955).
    EXPECT_LT(vaporOnlyDrift, 1e-6);
}

TEST(VaporTemperatureDiffusion, IsothermalControlAndGradient)
{
    // Control: uniform temperature -> the D1 term is inert, nothing should move.
    const auto isothermal = runClosedStrip([](double) { return 20.0; });

    // Test: linear temperature gradient 40 C (x=0) -> 20 C (x=L). Only the D1 vapour
    // temperature-gradient term can redistribute moisture.
    const auto gradient = runClosedStrip([](double x) { return 40.0 - 20.0 * (x / length); });

    std::cout << "[D1] isothermal drift = " << isothermal.maxDrift
              << ", gradient drift = " << gradient.maxDrift << "\n";

    // Isothermal: the D1 term is inert, so nothing moves and mass is exactly conserved.
    EXPECT_LT(isothermal.maxDrift, 1e-9);
    for(std::size_t col = 0; col < nCols; ++col)
    {
        EXPECT_NEAR(isothermal.finalBottomRow[col], phi0, 1e-12) << "column " << col;
    }

    // Gradient: the vapour temperature-gradient term drives moisture hot -> cold. The
    // closed-strip drift is bounded by the nonlinear-iteration residual only (measured
    // 2.3e-8 over 24 steps; the transport operator's column sums vanish).
    EXPECT_LT(gradient.maxDrift, 1e-6);

    // Day-end humidity checkpoints from the independently implemented 1D reference
    // solver (hygrothermfem_python, case vapor_gradient_cottaer) on the same grid,
    // steps and material tables. Tolerance ~2x the measured engine-reference
    // deviation of 4.1e-6.
    struct Checkpoint
    {
        std::size_t column;
        double humidity;
    };
    constexpr std::array<Checkpoint, 5> checkpoints{{{.column = 0, .humidity = 0.576382171496771},
                                                     {.column = 2, .humidity = 0.670102279072384},
                                                     {.column = 5, .humidity = 0.712475961763427},
                                                     {.column = 8, .humidity = 0.720351454778588},
                                                     {.column = 10, .humidity = 0.799275067054584}}};
    for(const auto & [column, humidity] : checkpoints)
    {
        EXPECT_NEAR(gradient.finalBottomRow[column], humidity, 1e-5) << "column " << column;
    }
}
