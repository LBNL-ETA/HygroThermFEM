#include <array>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

#include "Irregular1Mesh.hxx"
#include "TestMaterials.hxx"

// The THERM sample "Ireggular_1.thmz" (irregular three-material shape, quadtree mesh with
// skewed quadrilaterals) as a direct engine run. Documents the known-open start-up thermal
// ringing (doc/Moisture Governing Equations.md, "start-up thermal ringing"): ~48 of 143
// nodes flip-flop in temperature from step 2 (max ~9 K), damping to zero by step ~13, with
// the humidity field quiet. The flipping nodes cluster in Fiberglass Batts along the warm
// boundary; signature is step-count-anchored (same ringing at dt 360 and 36, smooth at
// 3.6) and mesh-refinement-independent, i.e. temporal-algorithmic. Suspected driver: the
// latent conductance coefficient k_lat = h_lg delta phi c_sat'(T) lagging the solve.
//
// DISABLED while the defect is parked. Promotion criterion once fixed: enable the test --
// the flip assertions at the bottom then guard the fix.
namespace
{
    struct FlipSummary
    {
        std::size_t totalFlips{0};
        double maxAmplitude{0.0};
        std::size_t worstStep{0};
        std::size_t lastFlipStep{0};
    };

    //! Temporal flip: a nodal value reversing direction between consecutive steps with
    //! both increments above the noise threshold. Physical relaxation under constant
    //! boundary conditions is monotone after the initial adjustment, so clustered large
    //! flips indicate solver flip-flop (same metric as D:/tmp/thmz_flip_detector.py).
    FlipSummary temporalFlips(const std::vector<std::vector<double>> & series,
                              const double threshold)
    {
        FlipSummary summary;
        if(series.size() < 3)
        {
            return summary;
        }
        const std::size_t nodeCount{series.front().size()};
        for(std::size_t node = 0; node < nodeCount; ++node)
        {
            for(std::size_t step = 0; step + 2 < series.size(); ++step)
            {
                const double diffFirst = series[step + 1][node] - series[step][node];
                const double diffSecond = series[step + 2][node] - series[step + 1][node];
                if(diffFirst * diffSecond >= 0.0)
                {
                    continue;
                }
                const double amplitude = (std::min)(std::abs(diffFirst), std::abs(diffSecond));
                if(amplitude <= threshold)
                {
                    continue;
                }
                ++summary.totalFlips;
                summary.lastFlipStep = (std::max)(summary.lastFlipStep, step + 1);
                if(amplitude > summary.maxAmplitude)
                {
                    summary.maxAmplitude = amplitude;
                    summary.worstStep = step + 1;
                }
            }
        }
        return summary;
    }

    void buildIrregular1Domain(HygroThermFEM::MultiDomain & multiDomain)
    {
        const HygroThermFEM::State initialState(
          {.temperature = 21.0, .humidity = 0.95, .pressure = 101325.0, .liquidPercent = 1.0});

        // Order matches Irregular1Mesh materialIndex: 0 = Fiberglass Batts,
        // 1 = Laminated panel, 2 = Stucco.
        const std::array<std::string, 3> materialNames{
          multiDomain.materials().createSolidMaterial(TestHelper::FiberglassBatts()).name(),
          multiDomain.materials().createSolidMaterial(TestHelper::LaminatedPanel()).name(),
          multiDomain.materials().createSolidMaterial(TestHelper::Stucco()).name()};

        for(const auto & node : Irregular1Mesh::nodes)
        {
            multiDomain.nodes().createNode(
              {.index = node.index, .x = node.xPos, .y = node.yPos, .state = initialState});
        }
        for(const auto & element : Irregular1Mesh::elements)
        {
            multiDomain.createElement({.node1 = element.nodeIDs[0],
                                       .node2 = element.nodeIDs[1],
                                       .node3 = element.nodeIDs[2],
                                       .node4 = element.nodeIDs[3],
                                       .material = materialNames[element.materialIndex]});
        }

        const HygroThermFEM::FixedBCHCCoefficients exterior{10.0, 7.0, 0.1};
        const HygroThermFEM::FixedBCHCCoefficients interior{40.0, 10.0, 0.3};
        for(const auto & edge : Irregular1Mesh::boundaryEdges)
        {
            multiDomain.createBC_FixedHc(edge.node1,
                                         edge.node2,
                                         edge.kind == Irregular1Mesh::EdgeKind::Indoor
                                           ? interior
                                           : exterior);
        }
    }
}   // namespace

TEST(ThermSample_Irregular1, DISABLED_StartupRinging)
{
    constexpr double dTime{360.0};
    constexpr unsigned numberOfSteps{100};

    HygroThermFEM::MultiDomain multiDomain;
    buildIrregular1Domain(multiDomain);

    auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
    auto humidities = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);
    std::vector<std::vector<double>> temperatureHistory{temperatures};
    std::vector<std::vector<double>> humidityHistory{humidities};

    for(unsigned step = 0; step < numberOfSteps; ++step)
    {
        const auto solution = multiDomain.transient(temperatures, humidities, dTime, step);
        temperatures = solution.temperature;
        humidities = solution.humidity;
        temperatureHistory.push_back(temperatures);
        humidityHistory.push_back(humidities);
    }

    const auto temperatureFlips = temporalFlips(temperatureHistory, 0.05);
    const auto humidityFlips = temporalFlips(humidityHistory, 0.002);

    std::cout << "[Irregular1] T flips: total=" << temperatureFlips.totalFlips
              << " maxAmp=" << temperatureFlips.maxAmplitude << " @step "
              << temperatureFlips.worstStep << " lastFlip@step "
              << temperatureFlips.lastFlipStep
              << " | phi flips: total=" << humidityFlips.totalFlips
              << " maxAmp=" << humidityFlips.maxAmplitude << std::endl;

    // Start-up ringing fix acceptance: no significant temperature direction reversals
    // beyond the first genuine adjustment step, and no large-amplitude flip anywhere.
    EXPECT_LT(temperatureFlips.maxAmplitude, 0.5);
    EXPECT_LE(temperatureFlips.lastFlipStep, 3u);
}
