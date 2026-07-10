#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "TestHelpers.hxx"
#include "TestMaterials.hxx"

using HygroThermFEM::State;

TEST(Moisture_2D_TwoElements_3, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Simple two elements example with moisture transfer.");

    HygroThermFEM::MultiDomain multiDomain({.performThermal = false, .performMoisture = true});

    const State state({
        .temperature = 20.0,
        .humidity = 0.0,
        .pressure = 101325.0,
        .liquidPercent = 1.0
    });

    multiDomain.nodes().createNode({.index = 1, .x = 0.15, .y = 0.05, .state = state});
    multiDomain.nodes().createNode({.index = 2, .x = 0.15, .y = 0, .state = state});
    multiDomain.nodes().createNode({.index = 3, .x = 0.05, .y = 0.05, .state = state});
    multiDomain.nodes().createNode({.index = 4, .x = 0.05, .y = 0, .state = state});
    multiDomain.nodes().createNode({.index = 5, .x = 0, .y = 0.05, .state = state});
    multiDomain.nodes().createNode({.index = 6, .x = 0, .y = 0, .state = state});

    const auto & material = multiDomain.materials().createSolidMaterial(TestHelper::CottaerSandstone());

    /// Create elements
    multiDomain.createElement({.node1 = 3, .node2 = 4, .node3 = 2, .node4 = 1, .material = material.name()});
    multiDomain.createElement({.node1 = 6, .node2 = 4, .node3 = 3, .node4 = 5, .material = material.name()});

    // Create Boundary Conditions
    constexpr auto airTemperature = 20.0;
    constexpr auto airHumidity = 0.5;

    const HygroThermFEM::TemperatureAndHumidity bcCoeff{airTemperature, airHumidity};

    multiDomain.moisture().createBC_FixedHumidity(5, 6, bcCoeff);

    constexpr auto dTime = 3600;
    constexpr auto nSteps = 24;

    auto humidities = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);
    std::vector<double> timesteps;
    std::vector<std::vector<double>> humiditySolution;
    std::vector<std::vector<HygroThermFEM::NodeFlux>> fluxSolution;

    for(unsigned i = 0; i < nSteps; ++i)
    {
        auto solution = multiDomain.moisture().transient(humidities, dTime).value();
        humidities = solution.solution;
        timesteps.push_back(solution.dTime);
        auto humidityContent = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);
        humiditySolution.push_back(humidityContent);
        fluxSolution.push_back(multiDomain.moisture().flux());
    }

    std::vector<double> correctTimesteps{3600, 3600, 3600, 3600, 3600, 3600, 3600, 3600,
                                         3600, 3600, 3600, 3600, 3600, 3600, 3600, 3600,
                                         3600, 3600, 3600, 3600, 3600, 3600, 3600, 3600};

    std::vector<std::vector<double>> correctHumiditySolution{{2.52320542e-06, 2.52320542e-06, 0.00129571287, 0.00129571287, 0.5, 0.5},
 {7.5549259e-06, 7.5549259e-06, 0.00258640516, 0.00258640516, 0.5, 0.5},
 {1.50805254e-05, 1.50805254e-05, 0.00387209009, 0.00387209009, 0.5, 0.5},
 {2.50854406e-05, 2.50854406e-05, 0.00515279035, 0.00515279035, 0.5, 0.5},
 {3.75551803e-05, 3.75551803e-05, 0.00642852848, 0.00642852848, 0.5, 0.5},
 {5.24753256e-05, 5.24753256e-05, 0.00769932697, 0.00769932697, 0.5, 0.5},
 {6.9831529e-05, 6.9831529e-05, 0.00896520817, 0.00896520817, 0.5, 0.5},
 {8.96095142e-05, 8.96095142e-05, 0.0102261943, 0.0102261943, 0.5, 0.5},
 {0.000111795076, 0.000111795076, 0.0114823076, 0.0114823076, 0.5, 0.5},
 {0.00013637408, 0.00013637408, 0.01273357, 0.01273357, 0.5, 0.5},
 {0.000163332461, 0.000163332461, 0.0139800035, 0.0139800035, 0.5, 0.5},
 {0.000192656225, 0.000192656225, 0.01522163, 0.01522163, 0.5, 0.5},
 {0.000224331447, 0.000224331447, 0.0164584711, 0.0164584711, 0.5, 0.5},
 {0.000258344272, 0.000258344272, 0.0176905484, 0.0176905484, 0.5, 0.5},
 {0.000294680911, 0.000294680911, 0.0189178836, 0.0189178836, 0.5, 0.5},
 {0.000333327649, 0.000333327649, 0.0201404981, 0.0201404981, 0.5, 0.5},
 {0.000374270833, 0.000374270833, 0.0213584131, 0.0213584131, 0.5, 0.5},
 {0.000417496882, 0.000417496882, 0.0225716499, 0.0225716499, 0.5, 0.5},
 {0.000462992283, 0.000462992283, 0.0237802296, 0.0237802296, 0.5, 0.5},
 {0.000510743586, 0.000510743586, 0.0249841733, 0.0249841733, 0.5, 0.5},
 {0.000560737413, 0.000560737413, 0.0261835018, 0.0261835018, 0.5, 0.5},
 {0.00061296045, 0.00061296045, 0.0273782361, 0.0273782361, 0.5, 0.5},
 {0.000667399448, 0.000667399448, 0.0285683968, 0.0285683968, 0.5, 0.5},
 {0.000724041227, 0.000724041227, 0.0297540045, 0.0297540045, 0.5, 0.5}};


    EXPECT_EQ(correctTimesteps.size(), timesteps.size());

    for(size_t i = 0u; i < correctTimesteps.size(); ++i)
    {
        EXPECT_NEAR(correctTimesteps[i], timesteps[i], 1e-6);
    }

    TestHelper::dumpGolden("correctHumiditySolution", humiditySolution);
    TestHelper::expectNear(correctHumiditySolution, humiditySolution, 1e-6);
}
