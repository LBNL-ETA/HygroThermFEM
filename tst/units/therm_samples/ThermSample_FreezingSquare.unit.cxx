#include <algorithm>
#include <cmath>
#include <numbers>
#include <vector>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

#include "FreezingSquareMesh.hxx"
#include "TestMaterials.hxx"

// The THERM benchmark "Freezing_Square_1.thmz" (0.5 x 0.2 m stucco slab, Simmetrix
// mesh, cold -10 C film-1000 edge, warm 0 C film-5 edge, fusion on, moisture off) as
// a direct engine run with THERM's GUI solver settings (relaxation 1, tolerance 1e-5,
// 25 iterations) AND the moisture/temperature-dependent conductivity enabled -- the
// stressed configuration that exposed the false-convergence energy leak: the tabular
// conductivity slows the fixed-point contraction, the change-metric acceptance then
// stagnates at large-residual iterates, and each accepted residual is a per-step
// energy imbalance (~6 % of the extracted heat, front lagging the Stefan solution by
// the same amount). Locks in the residual-reduction convergence that fusion runs use
// instead (ThermalDomain::useResidualConvergence).
namespace
{
    //! Solves lam * exp(lam^2) * erf(lam) = Ste / sqrt(pi) by bisection.
    double stefanLambda(const double stefanNumber)
    {
        const double target = stefanNumber / std::sqrt(std::numbers::pi);
        double low = 1.0e-6;
        double high = 3.0;
        for(int iteration = 0; iteration < 200; ++iteration)
        {
            const double mid = 0.5 * (low + high);
            const double value = mid * std::exp(mid * mid) * std::erf(mid) - target;
            if(value > 0.0)
            {
                high = mid;
            }
            else
            {
                low = mid;
            }
        }
        return 0.5 * (low + high);
    }

    //! Mid-band 1D profile along x from the unstructured node set; front where the
    //! temperature crosses `level` going up, linear interpolation.
    double frontPosition(const std::vector<double> & temperatures,
                         const double level)
    {
        std::vector<std::pair<double, double>> profile;
        for(const auto & node : FreezingSquareMesh::nodes)
        {
            if(std::abs(node.yPos) < 0.02)
            {
                profile.emplace_back(node.xPos, temperatures[node.index - 1]);
            }
        }
        std::sort(profile.begin(), profile.end());
        for(std::size_t idx = 1; idx < profile.size(); ++idx)
        {
            const auto & [xPrev, tPrev] = profile[idx - 1];
            const auto & [xCurr, tCurr] = profile[idx];
            if(tPrev < level && tCurr >= level)
            {
                const double fraction = (level - tPrev) / (tCurr - tPrev);
                return (xPrev + fraction * (xCurr - xPrev));
            }
        }
        return -1.0;
    }

    void buildFreezingSquareDomain(HygroThermFEM::MultiDomain & multiDomain)
    {
        const HygroThermFEM::State initialState(
          {.temperature = 0.0, .humidity = 0.99, .pressure = 101325.0, .liquidPercent = 1.0});

        const auto & stucco =
          multiDomain.materials().createSolidMaterial(TestHelper::Stucco());

        for(const auto & node : FreezingSquareMesh::nodes)
        {
            multiDomain.nodes().createNode(
              {.index = node.index, .x = node.xPos, .y = node.yPos, .state = initialState});
        }
        for(const auto & element : FreezingSquareMesh::elements)
        {
            multiDomain.createElement({.node1 = element.nodeIDs[0],
                                       .node2 = element.nodeIDs[1],
                                       .node3 = element.nodeIDs[2],
                                       .node4 = element.nodeIDs[3],
                                       .material = stucco.name()});
        }

        // Per-timestep coefficient vectors: the same path THERM's BCBuilder uses for the
        // transient BC series (constant -10 C / h=1000 cold, 0 C / h=5 warm, 240 entries).
        constexpr std::size_t seriesLength{240};
        const std::vector<HygroThermFEM::FixedBCHCCoefficients> coldSeries(
          seriesLength, HygroThermFEM::FixedBCHCCoefficients{-10.0, 1000.0, 0.8});
        const std::vector<HygroThermFEM::FixedBCHCCoefficients> warmSeries(
          seriesLength, HygroThermFEM::FixedBCHCCoefficients{0.0, 5.0, 0.99});
        for(const auto & edge : FreezingSquareMesh::boundaryEdges)
        {
            multiDomain.createBC_FixedHc(edge.node1,
                                         edge.node2,
                                         edge.kind == FreezingSquareMesh::EdgeKind::Cold
                                           ? coldSeries
                                           : warmSeries);
        }
    }
}   // namespace

TEST(ThermSample_FreezingSquare, StefanFront)
{
    HygroThermFEM::SimulationProperties::Instance().setCalculationParameters(
      true, true, true, true, true);
    HygroThermFEM::SimulationProperties::Instance().setExcludeLatentHeatOfFusion(false);
    HygroThermFEM::SimulationProperties::Instance().setIterationParameters(1.0, 1e-5, 25);

    HygroThermFEM::MultiDomain multiDomain;
    multiDomain.performThermalSimulation(true);
    multiDomain.performMoistureSimulation(false);
    buildFreezingSquareDomain(multiDomain);

    constexpr double dTime{600.0};
    constexpr unsigned numberOfSteps{144};

    auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
    auto humidities = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);
    std::vector<std::vector<double>> history{temperatures};
    for(unsigned step = 0; step < numberOfSteps; ++step)
    {
        const auto solution = multiDomain.transient(temperatures, humidities, dTime, step);
        temperatures = solution.temperature;
        humidities = solution.humidity;
        history.push_back(temperatures);
    }
    HygroThermFEM::SimulationProperties::Instance().reset();

    // Frozen-region properties as in Freezing.StefanFrontPosition (w = 95 at phi 0.99).
    const double waterContent{95.0};
    const double frozenCapacity =
      (1800.0 + waterContent) * (850.0 + waterContent / 916.7 * 2108.0);
    const double latentVolumetric = 333550.0 * waterContent;
    const double stefanNumber = frozenCapacity * 10.0 / latentVolumetric;
    const double lam = stefanLambda(stefanNumber);
    const double diffusivity = 0.85 / frozenCapacity;

    const double xStart = -0.25;
    for(const unsigned step : {48u, 96u, 144u})
    {
        const double analytic =
          2.0 * lam * std::sqrt(diffusivity * static_cast<double>(step) * dTime);
        const double numeric = frontPosition(history[step], -0.05) - xStart;
        EXPECT_NEAR(numeric, analytic, 0.05 * analytic) << "step " << step;
    }
}
