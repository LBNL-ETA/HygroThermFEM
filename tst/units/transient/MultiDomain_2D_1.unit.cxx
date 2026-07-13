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
 {0.00301405212073, 0.00301405212073, 1.06191775315, 1.06191775315, 2.11315044157, 2.11315044157},
 {0.00452702242003, 0.00452702242003, 1.06250629544, 1.06250629544, 2.1104603867, 2.1104603867},
 {0.00604247872518, 0.00604247872518, 1.06294061188, 1.06294061188, 2.10807629752, 2.10807629752},
 {0.00755948349239, 0.00755948349239, 1.06326290821, 1.06326290821, 2.10591470009, 2.10591470009},
 {0.00907724445524, 0.00907724445524, 1.06350276883, 1.06350276883, 2.10391721789, 2.10391721789},
 {0.0105951091997, 0.0105951091997, 1.06368134519, 1.06368134519, 2.10204220042, 2.10204220042},
 {0.0121125517642, 0.0121125517642, 1.0638140017, 1.0638140017, 2.10025944484, 2.10025944484},
 {0.0136291553339, 0.0136291553339, 1.06391204125, 1.06391204125, 2.09854676217, 2.09854676217},
 {0.0151445941226, 0.0151445941226, 1.06398386257, 1.06398386257, 2.09688768074, 2.09688768074}};

    TestHelper::dumpGolden("correctWaterContentSolution", waterContentSolution);
    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{{22.2951112517, 22.2951112517, 30.0026019148, 30.0026019148, 37.7098323131, 37.7098323131},
 {24.064706439, 24.064706439, 30.004501092, 30.004501092, 35.9443720234, 35.9443720234},
 {25.4291151299, 25.4291151299, 30.0059422183, 30.0059422183, 34.5832539879, 34.5832539879},
 {26.4811194012, 26.4811194012, 30.0070627531, 30.0070627531, 33.5337971575, 33.5337971575},
 {27.2922560303, 27.2922560303, 30.0079458848, 30.0079458848, 32.724601339, 32.724601339},
 {27.9176796984, 27.9176796984, 30.0086463036, 30.0086463036, 32.1006427827, 32.1006427827},
 {28.3999141716, 28.3999141716, 30.009202838, 30.009202838, 31.6195088584, 31.6195088584},
 {28.7717450015, 28.7717450015, 30.009644742, 30.009644742, 31.2485029193, 31.2485029193},
 {29.0584498894, 29.0584498894, 30.0099949244, 30.0099949244, 30.9624160851, 30.9624160851},
 {29.2795182323, 29.2795182323, 30.0102717094, 30.0102717094, 30.7418112167, 30.7418112167}};

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
 {0.00301405212073, 0.00301405212073, 1.06191775315, 1.06191775315, 2.11315044157, 2.11315044157},
 {0.00452702242003, 0.00452702242003, 1.06250629544, 1.06250629544, 2.1104603867, 2.1104603867},
 {0.00604247872518, 0.00604247872518, 1.06294061188, 1.06294061188, 2.10807629752, 2.10807629752},
 {0.00755948349239, 0.00755948349239, 1.06326290821, 1.06326290821, 2.10591470009, 2.10591470009},
 {0.00907724445524, 0.00907724445524, 1.06350276883, 1.06350276883, 2.10391721789, 2.10391721789},
 {0.0105951091997, 0.0105951091997, 1.06368134519, 1.06368134519, 2.10204220042, 2.10204220042},
 {0.0121125517642, 0.0121125517642, 1.0638140017, 1.0638140017, 2.10025944484, 2.10025944484},
 {0.0136291553339, 0.0136291553339, 1.06391204125, 1.06391204125, 2.09854676217, 2.09854676217},
 {0.0151445941226, 0.0151445941226, 1.06398386257, 1.06398386257, 2.09688768074, 2.09688768074}};

    TestHelper::dumpGolden("correctWaterContentSolution", waterContentSolution);
    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{{22.2951112517, 22.2951112517, 30.0026019148, 30.0026019148, 37.7098323131, 37.7098323131},
 {24.064706439, 24.064706439, 30.004501092, 30.004501092, 35.9443720234, 35.9443720234},
 {25.4291151299, 25.4291151299, 30.0059422183, 30.0059422183, 34.5832539879, 34.5832539879},
 {26.4811194012, 26.4811194012, 30.0070627531, 30.0070627531, 33.5337971575, 33.5337971575},
 {27.2922560303, 27.2922560303, 30.0079458848, 30.0079458848, 32.724601339, 32.724601339},
 {27.9176796984, 27.9176796984, 30.0086463036, 30.0086463036, 32.1006427827, 32.1006427827},
 {28.3999141716, 28.3999141716, 30.009202838, 30.009202838, 31.6195088584, 31.6195088584},
 {28.7717450015, 28.7717450015, 30.009644742, 30.009644742, 31.2485029193, 31.2485029193},
 {29.0584498894, 29.0584498894, 30.0099949244, 30.0099949244, 30.9624160851, 30.9624160851},
 {29.2795182323, 29.2795182323, 30.0102717094, 30.0102717094, 30.7418112167, 30.7418112167}};

    TestHelper::dumpGolden("correctTemperatureSolution", temperatureSolution);
    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-6);
}