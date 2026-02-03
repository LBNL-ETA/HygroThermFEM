#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::NodePool;
using HygroThermFEM::MaterialPool;

class Topaz2D_FluxBC : public testing::Test
{
protected:
    void SetUp() override
    {}

    void TearDown() override
    {
        NodePool::Instance().clear();
        MaterialPool::Instance().clear();
    }
};

TEST_F(Topaz2D_FluxBC, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Two elementsCreator example with transient.");

    // Enter nodes. Arguments are: node number, x-coordinate, y-coordinate, initial temperature

    // same temperature in every node (humidity and pressure irrelevant for this example)
    auto state = HygroThermFEM::State(0, 0, 101325, 0);

    NodePool::Instance().createNode(1, 0.15, 0.05, state);
    NodePool::Instance().createNode(2, 0.15, 0, state);
    NodePool::Instance().createNode(3, 0.05, 0.05, state);
    NodePool::Instance().createNode(4, 0.05, 0, state);
    NodePool::Instance().createNode(5, 0, 0.05, state);
    NodePool::Instance().createNode(6, 0, 0, state);

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

    auto & material =
      MaterialPool::Instance().createSolidMaterial("Cottaer Sandstone - non porous",
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

    HygroThermFEM::MultiDomain multiDomain(true, false);
    auto & domain = multiDomain.thermal();

    multiDomain.createElement(1, 2, 4, 3, material.name());
    multiDomain.createElement(5, 3, 4, 6, material.name());

    // Create Boundary Conditions
    // Positive flux means outside flow.
    constexpr auto surfaceFlux = -12.0;

    domain.createBC_FixedFlux(5, 6, surfaceFlux);

    constexpr auto dTime = 3600;
    constexpr auto nSteps = 4;

    auto temperatures = NodePool::Instance().properties(HygroThermFEM::Variable::temperature);
    std::vector<std::vector<double>> solution;

    for(unsigned i = 0; i < nSteps; ++i)
    {
        temperatures = domain.transient(temperatures, dTime).solution;
        solution.push_back(temperatures);
    }

    std::vector<std::vector<double>> correctSolution{
      {0.068797095, 0.068797095, 0.161296275, 0.161296275, 0.370195609, 0.370195609},
      {0.184225668, 0.184225668, 0.339421878, 0.339421878, 0.596640275, 0.596640275},
      {0.324790684, 0.324790684, 0.513783385, 0.513783385, 0.784104345, 0.784104345},
      {0.478007410, 0.478007410, 0.684010609, 0.684010609, 0.958667844, 0.958667844}};

    EXPECT_EQ(solution.size(), correctSolution.size());

    for(auto i = 0u; i < correctSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctSolution[i][j], solution[i][j], 1e-6);
        }
    }
}
