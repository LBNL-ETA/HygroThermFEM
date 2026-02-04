#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

/////////////////////////////////////////////////////////////////////////////////////
/// Transient temperature boundary conditions vs Analytical solution
///
/// This is test against analytical solution obtained from Carslaw-Jeager: page 122
/// NOTE: Carslaw-Jeager equation works only for specific coefficients (as used in example).
/////////////////////////////////////////////////////////////////////////////////////

class Analytical_ConvectionBC_Transient : public testing::Test
{
protected:
    void SetUp() override
    {}

    void TearDown() override
    {}
};

TEST_F(Analytical_ConvectionBC_Transient, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Example.");

    HygroThermFEM::MultiDomain multiDomain(true, false);

    // Enter nodes. Arguments are: node number, x-coordinate, y-coordinate, initial temperature

    /// Create slab that is 10 cm long and have nodes at every 1 cm
    const std::vector<double> gridXCoordinates{
      0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1};

    constexpr auto initialTemperature = 20.0;
    constexpr auto initialHumidity = 0.0;
    constexpr auto initialPressure = 101325.0;

    const HygroThermFEM::State state(initialTemperature, initialHumidity, initialPressure, 0);

    size_t nodeIndex = 0;
    for(const auto val : gridXCoordinates)
    {
        ++nodeIndex;
        multiDomain.nodes().createNode(nodeIndex, val, 0.00, state);
        ++nodeIndex;
        multiDomain.nodes().createNode(nodeIndex, val, 0.05, state);
    }

    // Material Properties (using C++20 designated initializers)
    const std::string materialName{"Test Material"};

    multiDomain.materials().createSolidMaterial({
        .name = materialName,
        .thermalConductivityDry = 1.8,
        .density = 2050.0,
        .porosity = 0.22,
        .heatCapacity = 850.0,
        .diffusionResistanceFactor = 15.0,
        .thermalConductivityMoistureDependent = {{0.0, 1.8}, {180, 1.8}},
        .moistureDependentMeasurementTemperature = 0.0,
        .thermalConductivityTemperatureDependent = {{0.0, 1.8}, {1, 1.8}},
        .temperatureDependentMeasurementHumidity = 0.0,
        .liquidTransportCurve = {{0, 0}, {180, 7E-7}},
        .sorptionCurve = {{0, 0}, {1, 180}}
    });

    /// Create elements
    for(size_t i = 1u; i <= (multiDomain.nodes().maxIndex() - 2) / 2; ++i)
    {
        const auto index1 = 2u * i - 1u;
        const auto index2 = 2u * i;
        const auto index3 = 2u * i + 2u;
        const auto index4 = 2u * i + 1u;

        multiDomain.createElement(index1, index2, index3, index4, materialName);
    }

    // Create Boundary Conditions
    constexpr auto tSurface = 0.0;

    multiDomain.thermal().createBC_FixedTemperature(21, 22, tSurface);

    constexpr auto dTime = 36;
    constexpr auto nSteps = 1000;

    auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
    std::vector<std::vector<double>> solution;

    for(unsigned i = 0; i < nSteps; ++i)
    {
        temperatures = multiDomain.thermal().transient(temperatures, dTime).solution;
        solution.push_back(temperatures);
    }

    std::vector<std::vector<double>> analyticalSolution{{10.171, 7.195, 0.000},
                                                        {4.064, 2.874, 0.000},
                                                        {1.623, 1.148, 0.000},
                                                        {0.649, 0.459, 0.000},
                                                        {0.259, 0.183, 0.000},
                                                        {0.104, 0.073, 0.000},
                                                        {0.041, 0.029, 0.000},
                                                        {0.017, 0.012, 0.000},
                                                        {0.007, 0.005, 0.000},
                                                        {0.003, 0.002, 0.000}};

    EXPECT_EQ(solution.size(), analyticalSolution.size() * 100);

    for(auto i = 0u; i < analyticalSolution.size(); ++i)
    {
        for(auto j = 0u; j < analyticalSolution[i].size(); ++j)
        {
            EXPECT_NEAR(analyticalSolution[i][j], solution[100 * i + 99][j * 10], 0.05);
        }
    }
}
