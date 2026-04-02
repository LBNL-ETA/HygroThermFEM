#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "TestHelpers.hxx"
#include "TestMaterials.hxx"

using HygroThermFEM::Nodes;
using HygroThermFEM::State;

////////////////////////////////////////////////////////////////////////////////////////////////////
/// This is simple two elements multi-domain example without boundary conditions. Initial
/// temperature and moisture distribution is not same in every node. This case should prove
/// that domain will try to reach equilibrium
////////////////////////////////////////////////////////////////////////////////////////////////////

TEST(MultiDomain_2D_1, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Simple two elements example with moisture and heat transfer.");

    HygroThermFEM::MultiDomain multiDomain;

    // Enter nodes. Arguments are: node number, x-coordinate, y-coordinate

    std::vector<double> gridXCoordinates{0, 0.05, 0.1};

    auto tVal = 0.0;
    auto deltaT = 10.0;
    auto hVal = 0.0;
    auto deltaH = 0.1;
    for(auto val : gridXCoordinates)
    {
        constexpr double initialPressure = 101325.0;
        constexpr double initialMoistureContent = 0.0;
        constexpr double initialTemperature = 20.0;

        multiDomain.nodes().createNode(
          {.x = val,
           .y = 0.00,
           .state = State{
              .temperature = initialTemperature + tVal,
              .humidity = initialMoistureContent + hVal,
              .pressure = initialPressure,
              .liquidPercent = 0
          }});
        multiDomain.nodes().createNode(
          {.x = val,
           .y = 0.05,
           .state = State{
              .temperature = initialTemperature + tVal,
              .humidity = initialMoistureContent + hVal,
              .pressure = initialPressure,
              .liquidPercent = 0
          }});
        tVal += deltaT;
        hVal += deltaH;
    }

    const auto & material = multiDomain.materials().createSolidMaterial(TestHelper::CottaerSandstone());

    /// Create elements
    for(size_t i = 1; i <= (multiDomain.nodes().maxIndex() - 2) / 2; ++i)
    {
        const auto node1 = 2u * i + 1u;
        const auto node2 = 2u * i + 2u;
        const auto node3 = 2u * i;
        const auto node4 = 2u * i - 1u;
        multiDomain.createElement({.node1 = node1, .node2 = node2, .node3 = node3, .node4 = node4, .material = material.name()});
    }

    constexpr auto dTime = 360;
    constexpr auto nSteps = 10;

    auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
    auto humidities = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);
    std::vector<std::vector<double>> temperatureSolution;
    std::vector<std::vector<double>> waterContentSolution;

    for(auto i = 0; i < nSteps; ++i)
    {
        auto aSolution = multiDomain.transient(temperatures, humidities, dTime);
        temperatureSolution.push_back(aSolution.temperature);
        waterContentSolution.push_back(aSolution.waterContent);
        temperatures = aSolution.temperature;
        humidities = aSolution.humidity;
    }

    std::vector<std::vector<double>> correctWaterContentSolution{
      {0.000827, 0.000827, 1.060090, 1.060090, 2.116574, 2.116574},
      {0.001771, 0.001771, 1.060141, 1.060141, 2.113687, 2.113687},
      {0.002814, 0.002814, 1.060169, 1.060169, 2.111179, 2.111179},
      {0.003939, 0.003939, 1.060184, 1.060184, 2.108945, 2.108945},
      {0.005130, 0.005130, 1.060190, 1.060190, 2.106912, 2.106912},
      {0.006375, 0.006375, 1.060190, 1.060190, 2.105027, 2.105027},
      {0.007662, 0.007662, 1.060188, 1.060188, 2.103253, 2.103253},
      {0.008981, 0.008981, 1.060184, 1.060184, 2.101563, 2.101563},
      {0.010326, 0.010326, 1.060179, 1.060179, 2.099938, 2.099938},
      {0.011690, 0.011690, 1.060173, 1.060173, 2.098361, 2.098361}};

    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{
      {22.291780, 22.291780, 30.001010, 30.001010, 37.719030, 37.719030},
      {24.058674, 24.058674, 30.002614, 30.002614, 35.958808, 35.958808},
      {25.421014, 25.421014, 30.004364, 30.004364, 34.600576, 34.600576},
      {26.471502, 26.471502, 30.006034, 30.006034, 33.552607, 33.552607},
      {27.281570, 27.281570, 30.007521, 30.007521, 32.744078, 32.744078},
      {27.906271, 27.906271, 30.008791, 30.008791, 32.120311, 32.120311},
      {28.388041, 28.388041, 30.009848, 30.009848, 31.639105, 31.639105},
      {28.759593, 28.759593, 30.010712, 30.010712, 31.267888, 31.267888},
      {29.046149, 29.046149, 30.011407, 30.011407, 30.981528, 30.981528},
      {29.267159, 29.267159, 30.011962, 30.011962, 30.760632, 30.760632}};

    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-6);
}

TEST(MultiDomain_2D_1, TestExample_1_Repeat)
{
    SCOPED_TRACE("Begin Test: Repeatability test.");

    HygroThermFEM::MultiDomain multiDomain;

    // Enter nodes. Arguments are: node number, x-coordinate, y-coordinate

    std::vector<double> gridXCoordinates{0, 0.05, 0.1};

    constexpr double initialTemperature = 20.0;
    constexpr double initialMoistureContent = 0.0;
    constexpr double initialPressure = 101325.0;

    auto tVal = 0.0;
    auto deltaT = 10.0;
    auto hVal = 0.0;
    auto deltaH = 0.1;
    for(auto val : gridXCoordinates)
    {
        multiDomain.nodes().createNode(
                {.x = val,
                 .y = 0.00,
                 .state = State{
                    .temperature = initialTemperature + tVal,
                    .humidity = initialMoistureContent + hVal,
                    .pressure = initialPressure,
                    .liquidPercent = 0
                }});
        multiDomain.nodes().createNode(
                {.x = val,
                 .y = 0.05,
                 .state = State{
                    .temperature = initialTemperature + tVal,
                    .humidity = initialMoistureContent + hVal,
                    .pressure = initialPressure,
                    .liquidPercent = 0
                }});
        tVal += deltaT;
        hVal += deltaH;
    }

    const auto & material = multiDomain.materials().createSolidMaterial(TestHelper::CottaerSandstone());

    /// Create elements
    for(size_t i = 1; i <= (multiDomain.nodes().maxIndex() - 2) / 2; ++i)
    {
        const auto node1 = 2u * i + 1u;
        const auto node2 = 2u * i + 2u;
        const auto node3 = 2u * i;
        const auto node4 = 2u * i - 1u;
        multiDomain.createElement({.node1 = node1, .node2 = node2, .node3 = node3, .node4 = node4, .material = material.name()});
    }

    constexpr auto dTime = 360;
    constexpr auto nSteps = 10;

    auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
    auto humidities = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);
    std::vector<std::vector<double>> temperatureSolution;
    std::vector<std::vector<double>> waterContentSolution;

    for(auto i = 0; i < nSteps; ++i)
    {
        auto aSolution = multiDomain.transient(temperatures, humidities, dTime);
        temperatureSolution.push_back(aSolution.temperature);
        waterContentSolution.push_back(aSolution.waterContent);
        temperatures = aSolution.temperature;
        humidities = aSolution.humidity;
    }

    std::vector<std::vector<double>> correctWaterContentSolution{
            {0.000827, 0.000827, 1.060090, 1.060090, 2.116574, 2.116574},
            {0.001771, 0.001771, 1.060141, 1.060141, 2.113687, 2.113687},
            {0.002814, 0.002814, 1.060169, 1.060169, 2.111179, 2.111179},
            {0.003939, 0.003939, 1.060184, 1.060184, 2.108945, 2.108945},
            {0.005130, 0.005130, 1.060190, 1.060190, 2.106912, 2.106912},
            {0.006375, 0.006375, 1.060190, 1.060190, 2.105027, 2.105027},
            {0.007662, 0.007662, 1.060188, 1.060188, 2.103253, 2.103253},
            {0.008981, 0.008981, 1.060184, 1.060184, 2.101563, 2.101563},
            {0.010326, 0.010326, 1.060179, 1.060179, 2.099938, 2.099938},
            {0.011690, 0.011690, 1.060173, 1.060173, 2.098361, 2.098361}};

    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{
            {22.291780, 22.291780, 30.001010, 30.001010, 37.719030, 37.719030},
            {24.058674, 24.058674, 30.002614, 30.002614, 35.958808, 35.958808},
            {25.421014, 25.421014, 30.004364, 30.004364, 34.600576, 34.600576},
            {26.471502, 26.471502, 30.006034, 30.006034, 33.552607, 33.552607},
            {27.281570, 27.281570, 30.007521, 30.007521, 32.744078, 32.744078},
            {27.906271, 27.906271, 30.008791, 30.008791, 32.120311, 32.120311},
            {28.388041, 28.388041, 30.009848, 30.009848, 31.639105, 31.639105},
            {28.759593, 28.759593, 30.010712, 30.010712, 31.267888, 31.267888},
            {29.046149, 29.046149, 30.011407, 30.011407, 30.981528, 30.981528},
            {29.267159, 29.267159, 30.011962, 30.011962, 30.760632, 30.760632}};

    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-6);
}