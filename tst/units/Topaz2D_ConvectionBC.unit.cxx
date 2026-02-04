#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

/////////////////////////////////////////////////////////////////////////////////////
/// Transient heat transfer example on Sandstone specimen using data from database
///   Lumped mass matrix
///   Time-step 1 hour
///   Six nodes block at initial temperatures in nodes of 100 degrees
///   Initial temperature boundary conditions at nodes 5 and 6 are 12 degrees
///   Solution achieved with linear solver (no iterations required in this case
/////////////////////////////////////////////////////////////////////////////////////

class Topaz2D_ConvectionBC : public testing::Test
{
protected:
    void SetUp() override
    {}

    void TearDown() override
    {}
};

/// Transient convection BC. Example is tested in EXCEL and Topaz2D.

TEST_F(Topaz2D_ConvectionBC, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Two elementsCreator example with transient.");

    HygroThermFEM::MultiDomain multiDomain(true, false);

    // Enter nodes. Arguments are: node number, x-coordinate, y-coordinate, initial temperature

    // same temperature in every node (humidity and pressure irrelevant for this example)
    const auto state = HygroThermFEM::State(0, 0, 101325, 0);

    multiDomain.nodePool().createNode(1, 0.15, 0.05, state);
    multiDomain.nodePool().createNode(2, 0.15, 0, state);
    multiDomain.nodePool().createNode(3, 0.05, 0.05, state);
    multiDomain.nodePool().createNode(4, 0.05, 0, state);
    multiDomain.nodePool().createNode(5, 0, 0.05, state);
    multiDomain.nodePool().createNode(6, 0, 0, state);

    // Material Properties (Cottaer Sandstone - non porous)
    constexpr double thermalConductivityDry{1.8};
    constexpr double density{2050.0};
    constexpr double porosity{0.0};
    constexpr double specificHeatCapacityDry{850.0};
    constexpr double diffusionResistanceFactor{15.0};
    const std::vector<FenestrationCommon::point> thermalConductivityMoistureDependent = {
      {0.0, 1.8}, {180, 1.8}};
    constexpr double thermalConductivityMeasuredAtTemperature{0};
    const std::vector<FenestrationCommon::point> thermalConductivityTemperatureDependent = {
      {0.0, 1.8}, {1, 1.8}};
    constexpr double thermalConductivityMeasuredAtHumidity{0};
    const std::vector<FenestrationCommon::point> liquidTransportationCurve = {{0, 0},
                                                                              {27, 1E-8},
                                                                              {45, 1.1E-8},
                                                                              {90, 2E-8},
                                                                              {126, 3.5E-8},
                                                                              {144, 5E-8},
                                                                              {162, 1E-7},
                                                                              {171, 2E-7},
                                                                              {180, 7E-7}};

    const std::vector<FenestrationCommon::point> moistureStorageFunction = {{0, 0},
                                                                            {0.5, 5.3},
                                                                            {0.65, 8.4},
                                                                            {0.8, 12},
                                                                            {0.93, 17},
                                                                            {0.95, 25},
                                                                            {0.99, 63},
                                                                            {0.995, 83},
                                                                            {0.999, 120},
                                                                            {1, 180}};

    const auto & material =
      multiDomain.materials().createSolidMaterial("Cottaer Sandstone - non porous",
                                                  thermalConductivityDry,
                                                  density,
                                                  porosity,
                                                  specificHeatCapacityDry,
                                                  diffusionResistanceFactor,
                                                  thermalConductivityMoistureDependent,
                                                  thermalConductivityMeasuredAtTemperature,
                                                  thermalConductivityTemperatureDependent,
                                                  thermalConductivityMeasuredAtHumidity,
                                                  liquidTransportationCurve,
                                                  moistureStorageFunction);

    multiDomain.createElement(3, 4, 2, 1, material.name());
    multiDomain.createElement(6, 4, 3, 5, material.name());

    // Create Boundary Conditions
    constexpr auto tSurface = 20.0;
    constexpr auto hc = 20.0;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff{tSurface, hc};

    multiDomain.thermal().createBC_FixedHc(5, 6, bcCoeff);

    constexpr auto dTime = 3600;
    constexpr auto nSteps = 4;


    auto temperatures = multiDomain.nodePool().properties(HygroThermFEM::Variable::temperature);
    std::vector<std::vector<double>> solution;

    for(unsigned i = 0; i < nSteps; ++i)
    {
        temperatures = multiDomain.thermal().transient(temperatures, dTime).solution;
        solution.push_back(temperatures);
    }

    std::vector<std::vector<double>> correctSolution{
      {1.4182108, 1.4182108, 3.3250259, 3.3250259, 7.6313602, 7.6313602},
      {3.4666897, 3.4666897, 6.2209137, 6.2209137, 10.518214, 10.518214},
      {5.6122051, 5.6122051, 8.4968969, 8.4968969, 12.234324, 12.234324},
      {7.6189241, 7.6189241, 10.317001, 10.317001, 13.501417, 13.501417}};

    EXPECT_EQ(solution.size(), correctSolution.size());

    for(auto i = 0u; i < correctSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctSolution[i][j], solution[i][j], 1e-6);
        }
    }
}
