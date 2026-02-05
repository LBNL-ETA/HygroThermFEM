#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "SlabCreator.hxx"

/////////////////////////////////////////////////////////////////////////////////////
/// Transient temperature boundary conditions vs Analytical solution
///
/// This is tested against analytical solution obtained from Carslaw-Jeager: page 97
/////////////////////////////////////////////////////////////////////////////////////

TEST(Analytical_TemperatureBC_Transient, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Example.");

    HygroThermFEM::MultiDomain multiDomain({.performThermal = true, .performMoisture = false});

    const HygroThermFEM::State state({
        .temperature = 1.0,
        .humidity = 0.0,
        .pressure = 101325.0,
        .liquidPercent = 0
    });

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

    /// Create slab that is 10 cm long and have nodes at every 1 cm
    TestHelper::SlabBuilder(multiDomain)
        .gridXCoordinates({0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0})
        .height(0.05)
        .material(material.name())
        .state(state)
        .build();

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
