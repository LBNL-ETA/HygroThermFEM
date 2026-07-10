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

    std::vector<std::vector<double>> correctWaterContentSolution{{0.0015046446, 0.0015046446, 1.06111329, 1.06111329, 2.11626878, 2.11626878},
 {0.00301390578, 0.00301390578, 1.06191905, 1.06191905, 2.113148, 2.113148},
 {0.00452669924, 0.00452669924, 1.06250941, 1.06250941, 2.11045448, 2.11045448},
 {0.00604200294, 0.00604200294, 1.06294572, 1.06294572, 2.10806656, 2.10806656},
 {0.00755890185, 0.00755890185, 1.06327001, 1.06327001, 2.10590108, 2.10590108},
 {0.00907660973, 0.00907660973, 1.0635118, 1.0635118, 2.10389979, 2.10389979},
 {0.0105944722, 0.0105944722, 1.06369221, 1.06369221, 2.10202111, 2.10202111},
 {0.0121119581, 0.0121119581, 1.06382659, 1.06382659, 2.10023485, 2.10023485},
 {0.013628644, 0.013628644, 1.06392627, 1.06392627, 2.09851882, 2.09851882},
 {0.015144198, 0.015144198, 1.06399964, 1.06399964, 2.09685652, 2.09685652}};

    TestHelper::dumpGolden("correctWaterContentSolution", waterContentSolution);
    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{{22.2917622, 22.2917622, 30.0009907, 30.0009907, 37.7191235, 37.7191235},
 {24.0586368, 24.0586368, 30.0025912, 30.0025912, 35.958948, 35.958948},
 {25.420961, 25.420961, 30.0043453, 30.0043453, 34.6007352, 34.6007352},
 {26.4714376, 26.4714376, 30.0060212, 30.0060212, 33.5527725, 33.5527725},
 {27.2814981, 27.2814981, 30.0075145, 30.0075145, 32.7442412, 32.7442412},
 {27.9061945, 27.9061945, 30.008791, 30.008791, 32.1204689, 32.1204689},
 {28.3879616, 28.3879616, 30.0098531, 30.0098531, 31.6392552, 31.6392552},
 {28.7595129, 28.7595129, 30.0107203, 30.0107203, 31.2680313, 31.2680313},
 {29.0460701, 29.0460701, 30.011419, 30.011419, 30.9816645, 30.9816645},
 {29.2670807, 29.2670807, 30.0119764, 30.0119764, 30.7607621, 30.7607621}};

    TestHelper::dumpGolden("correctTemperatureSolution", temperatureSolution);
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

    std::vector<std::vector<double>> correctWaterContentSolution{{0.0015046446, 0.0015046446, 1.06111329, 1.06111329, 2.11626878, 2.11626878},
 {0.00301390578, 0.00301390578, 1.06191905, 1.06191905, 2.113148, 2.113148},
 {0.00452669924, 0.00452669924, 1.06250941, 1.06250941, 2.11045448, 2.11045448},
 {0.00604200294, 0.00604200294, 1.06294572, 1.06294572, 2.10806656, 2.10806656},
 {0.00755890185, 0.00755890185, 1.06327001, 1.06327001, 2.10590108, 2.10590108},
 {0.00907660973, 0.00907660973, 1.0635118, 1.0635118, 2.10389979, 2.10389979},
 {0.0105944722, 0.0105944722, 1.06369221, 1.06369221, 2.10202111, 2.10202111},
 {0.0121119581, 0.0121119581, 1.06382659, 1.06382659, 2.10023485, 2.10023485},
 {0.013628644, 0.013628644, 1.06392627, 1.06392627, 2.09851882, 2.09851882},
 {0.015144198, 0.015144198, 1.06399964, 1.06399964, 2.09685652, 2.09685652}};

    TestHelper::dumpGolden("correctWaterContentSolution", waterContentSolution);
    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{{22.2917622, 22.2917622, 30.0009907, 30.0009907, 37.7191235, 37.7191235},
 {24.0586368, 24.0586368, 30.0025912, 30.0025912, 35.958948, 35.958948},
 {25.420961, 25.420961, 30.0043453, 30.0043453, 34.6007352, 34.6007352},
 {26.4714376, 26.4714376, 30.0060212, 30.0060212, 33.5527725, 33.5527725},
 {27.2814981, 27.2814981, 30.0075145, 30.0075145, 32.7442412, 32.7442412},
 {27.9061945, 27.9061945, 30.008791, 30.008791, 32.1204689, 32.1204689},
 {28.3879616, 28.3879616, 30.0098531, 30.0098531, 31.6392552, 31.6392552},
 {28.7595129, 28.7595129, 30.0107203, 30.0107203, 31.2680313, 31.2680313},
 {29.0460701, 29.0460701, 30.011419, 30.011419, 30.9816645, 30.9816645},
 {29.2670807, 29.2670807, 30.0119764, 30.0119764, 30.7607621, 30.7607621}};

    TestHelper::dumpGolden("correctTemperatureSolution", temperatureSolution);
    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-6);
}