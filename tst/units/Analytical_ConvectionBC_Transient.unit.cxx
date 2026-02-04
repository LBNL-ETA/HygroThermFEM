#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

/////////////////////////////////////////////////////////////////////////////////////
/// Transient temperature boundary conditions vs Analytical solution
///
/// This is test against analytical solution obtained from Carslaw-Jeager: page 97
/////////////////////////////////////////////////////////////////////////////////////

class Analytical_TemperatureBC_Transient : public testing::Test
{
protected:
    void SetUp() override
    {}

    void TearDown() override
    {}
};

TEST_F(Analytical_TemperatureBC_Transient, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Example.");

    HygroThermFEM::MultiDomain multiDomain(true, false);

    /// Create slab that is 10 cm long and have nodes at every 1 cm
    std::vector<double> gridXCoordinates{0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

    constexpr auto initialTemperature = 1.0;
    constexpr auto initialHumidity = 0.0;
    constexpr auto initialPressure = 101325.0;

    HygroThermFEM::State state(initialTemperature, initialHumidity, initialPressure, 0);

    size_t nodeIndex = 0;
    for(auto val : gridXCoordinates)
    {
        ++nodeIndex;
        multiDomain.nodes().createNode(nodeIndex, val, 0.00, state);
        ++nodeIndex;
        multiDomain.nodes().createNode(nodeIndex, val, 0.05, state);
    }

    // Material Properties (using C++20 designated initializers)
    const auto & material = multiDomain.materials().createSolidMaterial({
        .name = "Test Material",
        .thermalConductivityDry = 1.0,
        .density = 1.0,
        .porosity = 0.0,
        .heatCapacity = 1.0,
        .diffusionResistanceFactor = 15.0,
        .thermalConductivityMoistureDependent = {{0.0, 1.0}, {180, 1.0}},
        .moistureDependentMeasurementTemperature = 0.0,
        .thermalConductivityTemperatureDependent = {{0.0, 1.0}, {1, 1.0}},
        .temperatureDependentMeasurementHumidity = 0.0,
        .liquidTransportCurve = {{0, 0}, {180, 7E-7}},
        .sorptionCurve = {{0, 0}, {1, 180}}
    });

    /// Create elements
    for(size_t idx = 1; idx <= (multiDomain.nodes().maxIndex() - 2) / 2; ++idx)
    {
        const auto node1 = 2u * idx - 1u;
        const auto node2 = 2u * idx + 1u;
        const auto node3 = 2u * idx + 2u;
        const auto node4 = 2u * idx;

        multiDomain.createElement(node1, node2, node3, node4, material.name());
    }

    // Create Boundary Conditions
    constexpr auto tAir = 0.0;
    constexpr auto hc = 1.0;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff{tAir, hc};

    multiDomain.thermal().createBC_FixedHc(21, 22, bcCoeff);

    constexpr auto dTime = 0.001;
    constexpr auto nSteps = 1000;

    auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
    std::vector<std::vector<double>> solution;

    for(unsigned i = 0; i < nSteps; ++i)
    {
        temperatures = multiDomain.thermal().transient(temperatures, dTime).solution;
        solution.push_back(temperatures);
    }

    std::vector<std::vector<double>> analyticalSolution = {{0.99311, 0.95051, 0.72358},
                                                           {0.95064, 0.87925, 0.64339},
                                                           {0.89180, 0.81526, 0.58885},
                                                           {0.83095, 0.75671, 0.54417},
                                                           {0.77253, 0.70260, 0.50452},
                                                           {0.71768, 0.65243, 0.46827},
                                                           {0.66656, 0.60587, 0.43478},
                                                           {0.61903, 0.56264, 0.40374},
                                                           {0.57487, 0.52250, 0.37493},
                                                           {0.53386, 0.48522, 0.34818}};

    EXPECT_EQ(solution.size(), analyticalSolution.size() * 100);

    for(auto i = 0u; i < analyticalSolution.size(); ++i)
    {
        for(auto j = 0u; j < analyticalSolution[i].size(); ++j)
        {
            EXPECT_NEAR(analyticalSolution[i][j], solution[100 * i + 99][j * 10], 0.002);
        }
    }
}
