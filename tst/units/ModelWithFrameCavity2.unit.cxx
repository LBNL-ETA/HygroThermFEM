#include <gtest/gtest.h>
#include <memory>
#include <set>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::State;
using HygroThermFEM::ElementsLinear2D;
using HygroThermFEM::ElementThermalLinear2D;


TEST(TestModelWithFrameCavity2, TestDoubleFrameCavity)
{
    SCOPED_TRACE("Begin Test: Model with two frame cavities.");

    // Create elements grid
    HygroThermFEM::MultiDomain multiDomain(true, false);

    std::vector<double> gridX{0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07};
    std::vector<double> gridY{0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35};

    constexpr auto initialTemperature = 20;
    constexpr auto initialHumidity = 0.0;
    constexpr auto initialPressure = 101325.0;

    const State state(initialTemperature, initialHumidity, initialPressure);
    size_t nodeIndex = 0;

    // Crating grid nodes
    for(auto yVal : gridY)
    {
        for(auto xVal : gridX)
        {
            ++nodeIndex;
            multiDomain.nodes().createNode(nodeIndex, xVal, yVal, state);
        }
    }

    const auto & solidMaterial = multiDomain.materials().createSolidMaterial({
        .name = "Material 1",
        .thermalConductivityDry = 1.8,
        .density = 2050.0,
        .porosity = 0.22,
        .heatCapacity = 850.0,
        .diffusionResistanceFactor = 15.0,
        .thermalConductivityMoistureDependent = {{0.0, 1.8}, {180, 1.8}},
        .moistureDependentMeasurementTemperature = 0,
        .thermalConductivityTemperatureDependent = {{0.0, 1.8}, {1, 1.8}},
        .temperatureDependentMeasurementHumidity = 0,
        .liquidTransportCurve = {{0, 0}, {180, 2e-6}},
        .sorptionCurve = {{0, 0}, {1, 180}}
    });

    auto & frameCavity1 =
      multiDomain.materials().createGas("Frame Cavity 1", HygroThermFEM::CavityStandard::ISO15099);

    auto & frameCavity2 =
      multiDomain.materials().createGas("Frame Cavity 2", HygroThermFEM::CavityStandard::CEN);

    // Elements that will contain frame cavity
    std::set<size_t> frameCavity1Element{10, 11, 17, 18, 23, 24, 25, 30, 31, 32};
    std::set<size_t> frameCavity2Element{20, 21, 27, 28, 34, 35, 41, 42, 47, 48, 49};
    size_t elementNumber{0u};
    for(auto ix = 1u; ix < gridX.size(); ++ix)
    {
        for(auto iy = 1u; iy < gridY.size(); ++iy)
        {
            {
                ++elementNumber;
                const auto node1 = ix * gridX.size() + iy - gridX.size();
                const auto node2 = ix * gridX.size() + iy - gridX.size() + 1u;
                const auto node3 = ix * gridX.size() + (iy + 1u);
                const auto node4 = ix * gridX.size() + iy;
                std::string materialName;
                if(frameCavity1Element.find(elementNumber) != frameCavity1Element.end())
                {
                    materialName = frameCavity1.name();
                }
                else if(frameCavity2Element.find(elementNumber) != frameCavity2Element.end())
                {
                    materialName = frameCavity2.name();
                }
                else
                {
                    materialName = solidMaterial.name();
                }
                multiDomain.createElement(node1, node2, node3, node4, materialName);
            }
        }
    }

    // Create Boundary Conditions
    constexpr auto tAir = 0.0;
    constexpr auto hc = 30.0;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff{tAir, hc};

    // Build boundary condition nodes on left edge
    std::vector<size_t> bcnodes;
    for(size_t i = 0u; i < gridY.size(); ++i)
    {
        bcnodes.push_back(i * gridX.size() + 1);
    }

    // Now build boundary condition on left edge of domain rectangle
    for(size_t i = 1u; i < bcnodes.size(); ++i)
    {
        multiDomain.thermal().createBC_FixedHc(bcnodes[i - 1u], bcnodes[i], bcCoeff);
    }

    // Now perform transient calculation in order to make frame cavity update over the simulation
    constexpr auto dTime = 360;
    constexpr auto nSteps = 10;

    auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
    std::vector<std::vector<double>> solution;

    for(unsigned i = 0; i < nSteps; ++i)
    {
        temperatures = multiDomain.thermal().transient(temperatures, dTime).solution;
        solution.push_back(temperatures);
    }

    const auto correctThermalConductivity1{0.163195};
    const auto thermalCond1 = frameCavity1.thermalConductivityMoistureAndTemperatureDependent().getCurve()[0].y;
    EXPECT_NEAR(correctThermalConductivity1, thermalCond1, 1e-6);

    const auto correctThermalConductivity2{0.025013};
    const auto thermalCond2 = frameCavity2.thermalConductivityMoistureAndTemperatureDependent().getCurve()[0].y;
    EXPECT_NEAR(correctThermalConductivity2, thermalCond2, 1e-6);
}
