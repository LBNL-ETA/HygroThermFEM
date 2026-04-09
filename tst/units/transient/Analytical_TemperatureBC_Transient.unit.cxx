#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "SlabCreator.hxx"
#include "TestMaterials.hxx"

/////////////////////////////////////////////////////////////////////////////////////
/// Transient temperature boundary conditions vs Analytical solution
///
/// This is test against analytical solution obtained from Carslaw-Jeager: page 122
/// NOTE: Carslaw-Jeager equation works only for specific coefficients (as used in example).
/////////////////////////////////////////////////////////////////////////////////////

TEST(Analytical_ConvectionBC_Transient, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Example.");

    HygroThermFEM::MultiDomain multiDomain({.performThermal = true, .performMoisture = false});

    constexpr HygroThermFEM::State state({
        .temperature = 20.0,
        .humidity = 0.0,
        .pressure = 101325.0,
        .liquidPercent = 0
    });

    auto params = TestHelper::CottaerSandstone();
    params.name = "Test Material";
    params.liquidTransportCurve = {{0, 0}, {180, 7E-7}};
    params.sorptionCurve = {{0, 0}, {1, 180}};
    const auto & material = multiDomain.materials().createSolidMaterial(params);

    /// Create slab that is 10 cm long and have nodes at every 1 cm
    TestHelper::SlabBuilder(multiDomain)
        .gridXCoordinates({0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1})
        .height(0.05)
        .material(material.name())
        .state(state)
        .startCorner(TestHelper::StartCorner::BottomLeft)
        .direction(TestHelper::Direction::Clockwise)
        .build();

    // Create Boundary Conditions
    constexpr auto tSurface = 0.0;

    multiDomain.thermal().createBC_FixedTemperature(21, 22, tSurface);

    constexpr auto dTime = 36;
    constexpr auto nSteps = 1000;

    const auto solution = multiDomain.thermal().transientMultiStep(
      HygroThermFEM::Variable::temperature, dTime, nSteps);

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
