#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "SlabCreator.hxx"
#include "TestHelpers.hxx"
#include "TestMaterials.hxx"

using HygroThermFEM::Nodes;

TEST(MultiDomain_2D_2, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Simple two elements example with moisture and heat transfer.");

    HygroThermFEM::MultiDomain multiDomain;

    constexpr HygroThermFEM::State state({
        .temperature = 0.0,
        .humidity = 0.0,
        .pressure = 101325.0,
        .liquidPercent = 1.0
    });

    const auto & material = multiDomain.materials().createSolidMaterial(TestHelper::CottaerSandstone());

    TestHelper::SlabBuilder(multiDomain)
        .gridXCoordinates({0, 0.05, 0.1})
        .height(0.05)
        .material(material.name())
        .state(state)
        .startCorner(TestHelper::StartCorner::BottomRight)
        .direction(TestHelper::Direction::CounterClockwise)
        .build();

    /// Create Boundary Conditions
    constexpr auto hc = 1.0;
    constexpr auto airTemperature = 20.0;
    constexpr auto humidity = 0.6;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff{airTemperature, hc, humidity};

    multiDomain.createBC_FixedHc(1, 2, bcCoeff);

    constexpr auto dTime = 3600;
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

    std::vector<std::vector<double>> correctWaterContentSolution{{1.09461711009, 1.09461711009, 0.00106212274278, 0.00106212274278, 2.06118209937e-06, 2.06118209937e-06},
 {2.1245441741, 2.12454417409, 0.00332465868892, 0.00332465868893, 8.8360631258e-06, 8.8360631258e-06},
 {3.08858783046, 3.08858783045, 0.00681667559954, 0.00681667559955, 2.34888123003e-05, 2.34888123003e-05},
 {3.98615841199, 3.98615841198, 0.0115548475909, 0.0115548475909, 4.96366552083e-05, 4.96366552083e-05},
 {4.81708117862, 4.81708117861, 0.0175471888362, 0.0175471888362, 9.13186468369e-05, 9.13186468369e-05},
 {5.59094932774, 5.59094932773, 0.0246171513933, 0.0246171513933, 0.00015251985688, 0.00015251985688},
 {6.3231335437, 6.32313354367, 0.0324769936224, 0.0324769936224, 0.000236795164701, 0.000236795164701},
 {7.01334504449, 7.01334504445, 0.0411326478528, 0.0411326478528, 0.000347951055491, 0.00034795105549},
 {7.66162161962, 7.66162161958, 0.0505833111281, 0.0505833111282, 0.000490004716022, 0.000490004716022},
 {8.26825852107, 8.26825852104, 0.0608227761059, 0.060822776106, 0.000667148165096, 0.000667148165095}};

    TestHelper::dumpGolden("correctWaterContentSolution", waterContentSolution);
    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{{1.28169826644, 1.28169826644, 0.666638893203, 0.666638893203, 0.498933131692, 0.498933131692},
 {2.10758306379, 2.10758306379, 1.39459336614, 1.39459336614, 1.16927400258, 1.16927400258},
 {2.80351943401, 2.80351943401, 2.09888411237, 2.09888411237, 1.86502623422, 1.86502623422},
 {3.44167129511, 3.44167129511, 2.7682213349, 2.7682213349, 2.54101210619, 2.54101210619},
 {4.0393874807, 4.0393874807, 3.40176188482, 3.40176188482, 3.18523539655, 3.18523539655},
 {4.60677464765, 4.60677464765, 4.00279850868, 4.00279850868, 3.79714261856, 3.79714261856},
 {5.15227940322, 5.15227940322, 4.57673927576, 4.57673927576, 4.38064196592, 4.38064196592},
 {5.67346559886, 5.67346559886, 5.12486788461, 5.12486788461, 4.93767748202, 4.93767748202},
 {6.17034436023, 6.17034436023, 5.64788950728, 5.64788950728, 5.46926608957, 5.46926608957},
 {6.64349591694, 6.64349591694, 6.14648186746, 6.14648186746, 5.97617091747, 5.97617091747}};

    TestHelper::dumpGolden("correctTemperatureSolution", temperatureSolution);
    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-6);
}

TEST(MultiDomain_2D_2, TestExample_1_Repeat)
{
    SCOPED_TRACE("Begin Test: Repeatability test.");

    HygroThermFEM::MultiDomain multiDomain;

    constexpr HygroThermFEM::State state({
        .temperature = 0.0,
        .humidity = 0.0,
        .pressure = 101325.0,
        .liquidPercent = 1.0
    });

    const auto & material = multiDomain.materials().createSolidMaterial(TestHelper::CottaerSandstone());

    TestHelper::SlabBuilder(multiDomain)
        .gridXCoordinates({0, 0.05, 0.1})
        .height(0.05)
        .material(material.name())
        .state(state)
        .startCorner(TestHelper::StartCorner::BottomRight)
        .direction(TestHelper::Direction::CounterClockwise)
        .build();

    /// Create Boundary Conditions
    constexpr auto hc = 1.0;
    constexpr auto airTemperature = 20.0;
    constexpr auto humidity = 0.6;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff{airTemperature, hc, humidity};

    multiDomain.createBC_FixedHc(1, 2, bcCoeff);

    constexpr auto dTime = 3600;
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

    std::vector<std::vector<double>> correctWaterContentSolution{{1.09461711009, 1.09461711009, 0.00106212274278, 0.00106212274278, 2.06118209937e-06, 2.06118209937e-06},
 {2.1245441741, 2.12454417409, 0.00332465868892, 0.00332465868893, 8.8360631258e-06, 8.8360631258e-06},
 {3.08858783046, 3.08858783045, 0.00681667559954, 0.00681667559955, 2.34888123003e-05, 2.34888123003e-05},
 {3.98615841199, 3.98615841198, 0.0115548475909, 0.0115548475909, 4.96366552083e-05, 4.96366552083e-05},
 {4.81708117862, 4.81708117861, 0.0175471888362, 0.0175471888362, 9.13186468369e-05, 9.13186468369e-05},
 {5.59094932774, 5.59094932773, 0.0246171513933, 0.0246171513933, 0.00015251985688, 0.00015251985688},
 {6.3231335437, 6.32313354367, 0.0324769936224, 0.0324769936224, 0.000236795164701, 0.000236795164701},
 {7.01334504449, 7.01334504445, 0.0411326478528, 0.0411326478528, 0.000347951055491, 0.00034795105549},
 {7.66162161962, 7.66162161958, 0.0505833111281, 0.0505833111282, 0.000490004716022, 0.000490004716022},
 {8.26825852107, 8.26825852104, 0.0608227761059, 0.060822776106, 0.000667148165096, 0.000667148165095}};

    TestHelper::dumpGolden("correctWaterContentSolution", waterContentSolution);
    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{{1.28169826644, 1.28169826644, 0.666638893203, 0.666638893203, 0.498933131692, 0.498933131692},
 {2.10758306379, 2.10758306379, 1.39459336614, 1.39459336614, 1.16927400258, 1.16927400258},
 {2.80351943401, 2.80351943401, 2.09888411237, 2.09888411237, 1.86502623422, 1.86502623422},
 {3.44167129511, 3.44167129511, 2.7682213349, 2.7682213349, 2.54101210619, 2.54101210619},
 {4.0393874807, 4.0393874807, 3.40176188482, 3.40176188482, 3.18523539655, 3.18523539655},
 {4.60677464765, 4.60677464765, 4.00279850868, 4.00279850868, 3.79714261856, 3.79714261856},
 {5.15227940322, 5.15227940322, 4.57673927576, 4.57673927576, 4.38064196592, 4.38064196592},
 {5.67346559886, 5.67346559886, 5.12486788461, 5.12486788461, 4.93767748202, 4.93767748202},
 {6.17034436023, 6.17034436023, 5.64788950728, 5.64788950728, 5.46926608957, 5.46926608957},
 {6.64349591694, 6.64349591694, 6.14648186746, 6.14648186746, 5.97617091747, 5.97617091747}};

    TestHelper::dumpGolden("correctTemperatureSolution", temperatureSolution);
    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-6);
}