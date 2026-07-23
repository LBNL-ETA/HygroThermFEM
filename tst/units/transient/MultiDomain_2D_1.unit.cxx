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

    std::vector<std::vector<double>> correctWaterContentSolution{{0.00149299273008, 0.00149299273008, 1.06110992729, 1.06110992729, 2.11628715269, 2.11628715269},
 {0.00299512180912, 0.00299512180912, 1.06191289868, 1.06191289868, 2.11317908083, 2.11317908083},
 {0.00450360509443, 0.00450360509443, 1.06250078713, 1.06250078713, 2.11049482065, 2.11049482065},
 {0.00601632653004, 0.00601632653004, 1.06293482599, 1.06293482599, 2.10811402149, 2.10811402149},
 {0.00753168211332, 0.00753168211332, 1.0632570126, 1.0632570126, 2.10595429268, 2.10595429268},
 {0.00904845993815, 0.00904845993815, 1.0634968377, 1.0634968377, 2.10395786466, 2.10395786466},
 {0.0105657470057, 0.0105657470057, 1.06367541081, 1.06367541081, 2.10208343138, 2.10208343138},
 {0.0120828572436, 0.0120828572436, 1.06380807771, 1.06380807771, 2.10030098733, 2.10030098733},
 {0.013599276255, 0.013599276255, 1.06390613312, 1.06390613312, 2.0985884575, 2.0985884575},
 {0.0151146191515, 0.0151146191515, 1.06397797227, 1.06397797227, 2.09692943631, 2.09692943631}};

    TestHelper::dumpGolden("correctWaterContentSolution", waterContentSolution);
    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{{22.2951113267, 22.2951113267, 30.0026019635, 30.0026019635, 37.7098322881, 37.7098322881},
 {24.0647066093, 24.0647066093, 30.0045011986, 30.0045011986, 35.9443719697, 35.9443719697},
 {25.4291153929, 25.4291153929, 30.0059423782, 30.0059423782, 34.5832539069, 34.5832539069},
 {26.4811197456, 26.4811197456, 30.0070629581, 30.0070629581, 33.5337970524, 33.5337970524},
 {27.2922564429, 27.2922564429, 30.0079461264, 30.0079461264, 32.724601213, 32.724601213},
 {27.9176801664, 27.9176801664, 30.0086465742, 30.0086465742, 32.100642639, 32.100642639},
 {28.399914684, 28.399914684, 30.0092031314, 30.0092031314, 31.6195087002, 31.6195087002},
 {28.7717455491, 28.7717455491, 30.0096450531, 30.0096450531, 31.2485027492, 31.2485027492},
 {29.0584504646, 29.0584504646, 30.0099952493, 30.0099952493, 30.9624159055, 30.9624159055},
 {29.2795188291, 29.2795188291, 30.010272045, 30.010272045, 30.7418110295, 30.7418110295}};

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

    std::vector<std::vector<double>> correctWaterContentSolution{{0.00149299273008, 0.00149299273008, 1.06110992729, 1.06110992729, 2.11628715269, 2.11628715269},
 {0.00299512180912, 0.00299512180912, 1.06191289868, 1.06191289868, 2.11317908083, 2.11317908083},
 {0.00450360509443, 0.00450360509443, 1.06250078713, 1.06250078713, 2.11049482065, 2.11049482065},
 {0.00601632653004, 0.00601632653004, 1.06293482599, 1.06293482599, 2.10811402149, 2.10811402149},
 {0.00753168211332, 0.00753168211332, 1.0632570126, 1.0632570126, 2.10595429268, 2.10595429268},
 {0.00904845993815, 0.00904845993815, 1.0634968377, 1.0634968377, 2.10395786466, 2.10395786466},
 {0.0105657470057, 0.0105657470057, 1.06367541081, 1.06367541081, 2.10208343138, 2.10208343138},
 {0.0120828572436, 0.0120828572436, 1.06380807771, 1.06380807771, 2.10030098733, 2.10030098733},
 {0.013599276255, 0.013599276255, 1.06390613312, 1.06390613312, 2.0985884575, 2.0985884575},
 {0.0151146191515, 0.0151146191515, 1.06397797227, 1.06397797227, 2.09692943631, 2.09692943631}};

    TestHelper::dumpGolden("correctWaterContentSolution", waterContentSolution);
    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{{22.2951113267, 22.2951113267, 30.0026019635, 30.0026019635, 37.7098322881, 37.7098322881},
 {24.0647066093, 24.0647066093, 30.0045011986, 30.0045011986, 35.9443719697, 35.9443719697},
 {25.4291153929, 25.4291153929, 30.0059423782, 30.0059423782, 34.5832539069, 34.5832539069},
 {26.4811197456, 26.4811197456, 30.0070629581, 30.0070629581, 33.5337970524, 33.5337970524},
 {27.2922564429, 27.2922564429, 30.0079461264, 30.0079461264, 32.724601213, 32.724601213},
 {27.9176801664, 27.9176801664, 30.0086465742, 30.0086465742, 32.100642639, 32.100642639},
 {28.399914684, 28.399914684, 30.0092031314, 30.0092031314, 31.6195087002, 31.6195087002},
 {28.7717455491, 28.7717455491, 30.0096450531, 30.0096450531, 31.2485027492, 31.2485027492},
 {29.0584504646, 29.0584504646, 30.0099952493, 30.0099952493, 30.9624159055, 30.9624159055},
 {29.2795188291, 29.2795188291, 30.010272045, 30.010272045, 30.7418110295, 30.7418110295}};

    TestHelper::dumpGolden("correctTemperatureSolution", temperatureSolution);
    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-6);
}