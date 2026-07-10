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

    std::vector<std::vector<double>> correctWaterContentSolution{{0.00150464460051, 0.00150464460051, 1.06111328874, 1.06111328874, 2.11626877793, 2.11626877793},
 {0.00301401169943, 0.00301401169943, 1.06191835862, 1.06191835862, 2.11314927106, 2.11314927106},
 {0.00452694905178, 0.00452694905178, 1.06250774143, 1.06250774143, 2.11045756809, 2.11045756809},
 {0.00604240029877, 0.00604240029877, 1.06294296329, 1.06294296329, 2.10807167312, 2.10807167312},
 {0.00755943431857, 0.00755943431857, 1.06326615359, 1.06326615359, 2.1059082585, 2.1059082585},
 {0.00907725848842, 0.00907725848842, 1.06350686515, 1.06350686515, 2.10390901121, 2.10390901121},
 {0.0105952170688, 0.0105952170688, 1.06368623928, 1.06368623928, 2.10203230438, 2.10203230438},
 {0.0121127797661, 0.0121127797661, 1.06381964028, 1.06381964028, 2.10024793968, 2.10024793968},
 {0.0136295253932, 0.0136295253932, 1.06391837529, 1.06391837529, 2.09853372402, 2.09853372402},
 {0.0151451241857, 0.0151451241857, 1.0639908489, 1.0639908489, 2.09687317801, 2.09687317801}};

    TestHelper::dumpGolden("correctWaterContentSolution", waterContentSolution);
    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{{22.2940051262, 22.2940051262, 30.0021573416, 30.0021573416, 37.7143991, 37.7143991},
 {24.0627091868, 24.0627091868, 30.0041546329, 30.0041546329, 35.9514983929, 35.9514983929},
 {25.4264531103, 25.4264531103, 30.0059222798, 30.0059222798, 34.5917660855, 34.5917660855},
 {26.477991442, 26.477991442, 30.0074379544, 30.0074379544, 33.5430090423, 33.5430090423},
 {27.288821609, 27.288821609, 30.0087079073, 30.0087079073, 32.7341156724, 32.7341156724},
 {27.9140594531, 27.9140594531, 30.0097536914, 30.0097536914, 32.1102339725, 32.1102339725},
 {28.3961953077, 28.3961953077, 30.0106035079, 30.0106035079, 31.6290536754, 31.6290536754},
 {28.7679879889, 28.7679879889, 30.0112869636, 30.0112869636, 31.257938901, 31.257938901},
 {29.0546948266, 29.0546948266, 30.0118321508, 30.0118321508, 30.9717164326, 30.9717164326},
 {29.2757902572, 29.2757902572, 30.0122642185, 30.0122642185, 30.7509695913, 30.7509695913}};

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

    std::vector<std::vector<double>> correctWaterContentSolution{{0.00150464460051, 0.00150464460051, 1.06111328874, 1.06111328874, 2.11626877793, 2.11626877793},
 {0.00301401169943, 0.00301401169943, 1.06191835862, 1.06191835862, 2.11314927106, 2.11314927106},
 {0.00452694905178, 0.00452694905178, 1.06250774143, 1.06250774143, 2.11045756809, 2.11045756809},
 {0.00604240029877, 0.00604240029877, 1.06294296329, 1.06294296329, 2.10807167312, 2.10807167312},
 {0.00755943431857, 0.00755943431857, 1.06326615359, 1.06326615359, 2.1059082585, 2.1059082585},
 {0.00907725848842, 0.00907725848842, 1.06350686515, 1.06350686515, 2.10390901121, 2.10390901121},
 {0.0105952170688, 0.0105952170688, 1.06368623928, 1.06368623928, 2.10203230438, 2.10203230438},
 {0.0121127797661, 0.0121127797661, 1.06381964028, 1.06381964028, 2.10024793968, 2.10024793968},
 {0.0136295253932, 0.0136295253932, 1.06391837529, 1.06391837529, 2.09853372402, 2.09853372402},
 {0.0151451241857, 0.0151451241857, 1.0639908489, 1.0639908489, 2.09687317801, 2.09687317801}};

    TestHelper::dumpGolden("correctWaterContentSolution", waterContentSolution);
    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{{22.2940051262, 22.2940051262, 30.0021573416, 30.0021573416, 37.7143991, 37.7143991},
 {24.0627091868, 24.0627091868, 30.0041546329, 30.0041546329, 35.9514983929, 35.9514983929},
 {25.4264531103, 25.4264531103, 30.0059222798, 30.0059222798, 34.5917660855, 34.5917660855},
 {26.477991442, 26.477991442, 30.0074379544, 30.0074379544, 33.5430090423, 33.5430090423},
 {27.288821609, 27.288821609, 30.0087079073, 30.0087079073, 32.7341156724, 32.7341156724},
 {27.9140594531, 27.9140594531, 30.0097536914, 30.0097536914, 32.1102339725, 32.1102339725},
 {28.3961953077, 28.3961953077, 30.0106035079, 30.0106035079, 31.6290536754, 31.6290536754},
 {28.7679879889, 28.7679879889, 30.0112869636, 30.0112869636, 31.257938901, 31.257938901},
 {29.0546948266, 29.0546948266, 30.0118321508, 30.0118321508, 30.9717164326, 30.9717164326},
 {29.2757902572, 29.2757902572, 30.0122642185, 30.0122642185, 30.7509695913, 30.7509695913}};

    TestHelper::dumpGolden("correctTemperatureSolution", temperatureSolution);
    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-6);
}