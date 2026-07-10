#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "TestMaterials.hxx"
#include "TestHelpers.hxx"
#include "ObserveSimulationProgress.hxx"

using HygroThermFEM::Nodes;

TEST(MultiDomain_HighHumidity, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Simple two elements example with moisture transfer.");

    HygroThermFEM::MultiDomain multiDomain;

    std::vector<double> gridXCoordinates{0.15, 0.05, 0.00};

    constexpr HygroThermFEM::State state(
      {.temperature = 0.0, .humidity = 0.999, .pressure = 101325.0, .liquidPercent = 1.0});

    for(auto val : gridXCoordinates)
    {
        multiDomain.nodes().createNode({.x = val, .y = 0.05, .state = state});
        multiDomain.nodes().createNode({.x = val, .y = 0.00, .state = state});
    }

    const auto & material =
      multiDomain.materials().createSolidMaterial(TestHelper::CottaerSandstone());

    TestHelper::ObserveSimulationProgress progressThermal;
    multiDomain.subscribeThermal(&progressThermal);

    TestHelper::ObserveSimulationProgress progressMoisture;
    multiDomain.subscribeMoisture(&progressMoisture);

    /// Create elements
    for(size_t i = 1; i <= (multiDomain.nodes().maxIndex() - 2) / 2; ++i)
    {
        const auto node1 = 2u * i + 1u;
        const auto node2 = 2u * i + 2u;
        const auto node3 = 2u * i;
        const auto node4 = 2u * i - 1u;
        multiDomain.createElement({.node1 = node2,
                                   .node2 = node3,
                                   .node3 = node4,
                                   .node4 = node1,
                                   .material = material.name()});
    }

    // Create Boundary Conditions
    constexpr auto hc = 10.0;
    constexpr auto airTemperature = 20.0;
    constexpr auto airHumidity = 1.0;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff{airTemperature, hc, airHumidity};

    multiDomain.createBC_FixedHc(5, 6, bcCoeff);

    const auto results = multiDomain.transientMultiStep(3600, 10);

    const std::vector<double> correctHumidityError{5.17474753744e-11, 3.83729569973e-11, 2.94630418403e-11, 2.23246732669e-11, 1.68449918764e-11, 1.26434745926e-11, 9.45327600077e-12, 7.05477994472e-12, 5.25652644276e-12, 3.90172983404e-12};
    const std::vector<std::vector<double>> correctWaterContentSolution{{120.006901846, 120.006901846, 120.367262911, 120.367262911, 132.706612558, 132.706612558},
 {120.027040467, 120.027040467, 121.033181415, 121.033181415, 140.720396869, 140.720396869},
 {120.065158984, 120.065158984, 121.950508047, 121.950508047, 146.144529307, 146.144529307},
 {120.125336604, 120.125336604, 123.080088456, 123.080088455, 149.805192711, 149.805192711},
 {120.210161907, 120.210161907, 124.33538634, 124.33538634, 152.202534498, 152.202534498},
 {120.320963395, 120.320963395, 125.645946908, 125.645946908, 153.673178417, 153.673178417},
 {120.458290463, 120.458290463, 126.956591219, 126.956591219, 154.471762231, 154.471762231},
 {120.621758106, 120.621758106, 128.226640266, 128.226640266, 154.793078978, 154.793078978},
 {120.810208039, 120.810208039, 129.429000126, 129.429000126, 154.781301902, 154.781301902},
 {121.021990329, 121.021990329, 130.547746926, 130.547746926, 154.539981634, 154.539981634}};

    TestHelper::dumpGolden("correctHumidityError", results.moisture.errors);
    TestHelper::expectNear(correctHumidityError, results.moisture.errors, 1e-6);
    TestHelper::dumpGolden("correctWaterContentSolution", results.moisture.values);
    TestHelper::expectNear(correctWaterContentSolution, results.moisture.values, 1e-6);

    const std::vector<double> correctTemperatureError{6.44252504273e-07, 2.04956097033e-07, 1.30007396018e-07, 1.00712370583e-07, 8.22748020877e-08, 6.8613744217e-08, 5.78684707596e-08, 4.91719029449e-08, 4.20158205397e-08, 3.60621524647e-08};
    const std::vector<std::vector<double>> correctTemperatureSolution{{0.798766489892, 0.798766489892, 2.60750150089, 2.60750150089, 7.95260916693, 7.95260916693},
 {2.1084578976, 2.1084578976, 5.07398471509, 5.07398471509, 10.7644618285, 10.7644618285},
 {3.63443059712, 3.63443059712, 7.09086476445, 7.09086476445, 12.2674515859, 12.2674515859},
 {5.19881366865, 5.19881366865, 8.74437841012, 8.74437841012, 13.3508552162, 13.3508552162},
 {6.70768282238, 6.70768282238, 10.1300854975, 10.1300854975, 14.2219477446, 14.2219477446},
 {8.11609437079, 8.11609437079, 11.3135123466, 11.3135123466, 14.9506672609, 14.9506672609},
 {9.40610762087, 9.40610762087, 12.33760544, 12.33760544, 15.5721452658, 15.5721452658},
 {10.5742344098, 10.5742344098, 13.2314873455, 13.2314873455, 16.1083017668, 16.1083017668},
 {11.6244449404, 11.6244449404, 14.0159911291, 14.0159911291, 16.5742859166, 16.5742859166},
 {12.564305752, 12.564305752, 14.7068503662, 14.7068503662, 16.9812577733, 16.9812577733}};

    TestHelper::dumpGolden("correctTemperatureError", results.temperature.errors);
    TestHelper::expectNear(correctTemperatureError, results.temperature.errors, 1e-6);
    TestHelper::dumpGolden("correctTemperatureSolution", results.temperature.values);
    TestHelper::expectNear(correctTemperatureSolution, results.temperature.values, 1e-6);

    // Checking number of iterations within subiterations
    EXPECT_EQ(progressMoisture.calls(), (std::vector<unsigned>{}));
    EXPECT_EQ(progressThermal.calls(), (std::vector<unsigned>{}));
}

TEST(MultiDomain_HighHumidity, HighHumidityAndTemperature)
{
    SCOPED_TRACE("Begin Test: Simple two elements example with moisture transfer.");

    HygroThermFEM::MultiDomain multiDomain;

    const std::vector gridXCoordinates{0.15, 0.05, 0.00};

    constexpr HygroThermFEM::State state({
        .temperature = 30.0,
        .humidity = 0.9999,
        .pressure = 101325.0,
        .liquidPercent = 1.0
    });

    for(auto val : gridXCoordinates)
    {
        multiDomain.nodes().createNode({.x = val, .y = 0.05, .state = state});
        multiDomain.nodes().createNode({.x = val, .y = 0.00, .state = state});
    }

    const auto & material =
      multiDomain.materials().createSolidMaterial(TestHelper::CottaerSandstone());

    TestHelper::ObserveSimulationProgress progressThermal;
    multiDomain.subscribeThermal(&progressThermal);

    TestHelper::ObserveSimulationProgress progressMoisture;
    multiDomain.subscribeMoisture(&progressMoisture);

    /// Create elements
    for(size_t i = 1; i <= (multiDomain.nodes().maxIndex() - 2) / 2; ++i)
    {
        const auto node1 = 2u * i + 1u;
        const auto node2 = 2u * i + 2u;
        const auto node3 = 2u * i;
        const auto node4 = 2u * i - 1u;
        multiDomain.createElement({.node1 = node2,
                                   .node2 = node3,
                                   .node3 = node4,
                                   .node4 = node1,
                                   .material = material.name()});
    }

    // Create Boundary Conditions
    constexpr auto hc = 10.0;
    constexpr auto airTemperature = 20.0;
    constexpr auto airHumidity = 1.0;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff{airTemperature, hc, airHumidity};

    multiDomain.createBC_FixedHc(5, 6, bcCoeff);

    const auto results = multiDomain.transientMultiStep(3600, 10);

    const std::vector<double> correctHumidityError{1.23486463579e-11, 1.20713899358e-11, 9.84114983682e-12, 8.01520868319e-12, 6.68118408776e-12, 5.64869383174e-12, 4.81710703127e-12, 4.13311707688e-12, 3.55882846431e-12, 3.06885652104e-12};
    const std::vector<std::vector<double>> correctWaterContentSolution{{173.775711862, 173.775711862, 172.633351641, 172.633351641, 164.054650143, 164.054650143},
 {173.431967432, 173.431967432, 171.518921463, 171.518921463, 162.888035565, 162.888035565},
 {173.022354694, 173.022354694, 170.560097895, 170.560097895, 161.843817361, 161.843817361},
 {172.575593277, 172.575593277, 169.713402994, 169.713402994, 161.181812328, 161.181812328},
 {172.113909398, 172.113909398, 168.967202102, 168.967202102, 160.763150001, 160.763150001},
 {171.651613742, 171.651613742, 168.308815918, 168.308815918, 160.507945153, 160.507945153},
 {171.197688384, 171.197688384, 167.726832334, 167.726832334, 160.367269493, 160.367269493},
 {170.754509457, 170.754509457, 167.213412467, 167.213412467, 160.310009389, 160.310009389},
 {170.323727971, 170.323727971, 166.76124877, 166.76124877, 160.315033714, 160.315033714},
 {169.908409716, 169.908409716, 166.362593562, 166.362593562, 160.366733076, 160.366733076}};

    TestHelper::dumpGolden("correctHumidityError", results.moisture.errors);
    TestHelper::expectNear(correctHumidityError, results.moisture.errors, 1e-6);
    TestHelper::dumpGolden("correctWaterContentSolution", results.moisture.values);
    TestHelper::expectNear(correctWaterContentSolution, results.moisture.values, 1e-6);

    const std::vector<double> correctTemperatureError{1.88846121105e-08, 2.39467418397e-08, 2.22009965927e-08, 2.04867618903e-08, 1.88269709367e-08, 1.72451179876e-08, 1.57555729228e-08, 1.43637051623e-08, 1.30715694433e-08, 1.18785362079e-08};
    const std::vector<std::vector<double>> correctTemperatureSolution{{29.5373563673, 29.5373563673, 28.2848499156, 28.2848499156, 24.166636212, 24.166636212},
 {28.8796732398, 28.8796732398, 27.1053442647, 27.1053442647, 23.8330526625, 23.8330526625},
 {28.1508571465, 28.1508571465, 26.1881047177, 26.1881047177, 23.3535658046, 23.3535658046},
 {27.417501891, 27.417501891, 25.4457267327, 25.4457267327, 22.9620024098, 22.9620024098},
 {26.7140711146, 26.7140711146, 24.8256789631, 24.8256789631, 22.6325128729, 22.6325128729},
 {26.0569750011, 26.0569750011, 24.2956298768, 24.2956298768, 22.3492138443, 22.3492138443},
 {25.4526185207, 25.4526185207, 23.835024356, 23.835024356, 22.1018678083, 22.1018678083},
 {24.9020239194, 24.9020239194, 23.4302916154, 23.4302916154, 21.8836571445, 21.8836571445},
 {24.4033477413, 24.4033477413, 23.0720113204, 23.0720113204, 21.6898204521, 21.6898204521},
 {23.9533242655, 23.9533242655, 22.7532883374, 22.7532883374, 21.5168554009, 21.5168554009}};

    TestHelper::dumpGolden("correctTemperatureError", results.temperature.errors);
    TestHelper::expectNear(correctTemperatureError, results.temperature.errors, 1e-6);
    TestHelper::dumpGolden("correctTemperatureSolution", results.temperature.values);
    TestHelper::expectNear(correctTemperatureSolution, results.temperature.values, 1e-6);

    // Checking number of iterations within subiterations
    EXPECT_EQ(progressMoisture.calls(), (std::vector<unsigned>{}));
    EXPECT_EQ(progressThermal.calls(), (std::vector<unsigned>{}));
}

TEST(MultiDomain_HighHumidity, ExtremeHumidityAndTemperature)
{
    SCOPED_TRACE("Begin Test: Simple two elements example with moisture transfer.");

    HygroThermFEM::MultiDomain multiDomain;

    std::vector<double> gridXCoordinates{0.15, 0.05, 0.00};

    constexpr HygroThermFEM::State state({
        .temperature = 80.0,
        .humidity = 0.9999,
        .pressure = 101325.0,
        .liquidPercent = 1.0
    });

    for(auto val : gridXCoordinates)
    {
        multiDomain.nodes().createNode({.x = val, .y = 0.05, .state = state});
        multiDomain.nodes().createNode({.x = val, .y = 0.00, .state = state});
    }

    const auto & material =
      multiDomain.materials().createSolidMaterial(TestHelper::CottaerSandstone());

    TestHelper::ObserveSimulationProgress progressThermal;
    multiDomain.subscribeThermal(&progressThermal);

    TestHelper::ObserveSimulationProgress progressMoisture;
    multiDomain.subscribeMoisture(&progressMoisture);

    /// Create elements
    for(size_t i = 1; i <= (multiDomain.nodes().maxIndex() - 2) / 2; ++i)
    {
        const auto node1 = 2u * i + 1u;
        const auto node2 = 2u * i + 2u;
        const auto node3 = 2u * i;
        const auto node4 = 2u * i - 1u;
        multiDomain.createElement({.node1 = node2,
                                   .node2 = node3,
                                   .node3 = node4,
                                   .node4 = node1,
                                   .material = material.name()});
    }

    // Create Boundary Conditions
    constexpr auto hc = 10.0;
    constexpr auto airTemperature = 20.0;
    constexpr auto airHumidity = 1.0;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff{airTemperature, hc, airHumidity};

    multiDomain.createBC_FixedHc(5, 6, bcCoeff);

    const auto results = multiDomain.transientMultiStep(3600, 10);

    const std::vector<double> correctHumidityError{1.68966287155e-08, 3.97584426105e-09, 3.21510322436e-10, 6.98879334584e-11, 1.42544501188e-10, 9.0492231366e-11, 5.94894845825e-11, 2.11681241051e-11, 2.22305269707e-10, 2.71388114025e-10};
    const std::vector<std::vector<double>> correctWaterContentSolution{{170.641312207, 170.641312207, 129.709323445, 129.709323445, 25.7909135045, 25.7909135045},
 {167.602024664, 167.602024664, 124.4237794, 124.423779399, 66.1425674603, 66.1425674603},
 {165.031880127, 165.031880127, 120.754819413, 120.754819413, 72.6750421244, 72.6750421244},
 {162.312054612, 162.312054612, 117.499838997, 117.499838997, 68.3934610936, 68.3934610936},
 {159.42511222, 159.42511222, 115.042150509, 115.042150509, 66.8468616519, 66.8468616519},
 {156.521102508, 156.521102508, 113.064076515, 113.064076515, 65.0318129675, 65.0318129675},
 {153.693197413, 153.693197413, 111.446461902, 111.446461902, 63.7938384997, 63.7938384997},
 {150.985859039, 150.985859039, 110.09405251, 110.09405251, 62.9581053047, 62.9581053047},
 {148.42949563, 148.42949563, 108.90041801, 108.90041801, 62.645516771, 62.645516771},
 {146.037392681, 146.037392681, 107.887541879, 107.887541879, 62.5844329003, 62.5844329003}};

    TestHelper::dumpGolden("correctHumidityError", results.moisture.errors);
    TestHelper::expectNear(correctHumidityError, results.moisture.errors, 1e-6);
    TestHelper::dumpGolden("correctWaterContentSolution", results.moisture.values);
    TestHelper::expectNear(correctWaterContentSolution, results.moisture.values, 1e-6);

    const std::vector<double> correctTemperatureError{9.85693153543e-07, 2.33324258805e-09, 8.46988647312e-08, 5.29755991827e-08, 5.60510345442e-08, 5.02672645921e-08, 4.7065963295e-08, 4.35336315818e-08, 4.02553492046e-08, 3.71014339597e-08};
    const std::vector<std::vector<double>> correctTemperatureSolution{{71.4425593867, 71.4425593867, 48.594921971, 48.594921971, -24.7970870229, -24.7970870229},
 {63.9418923073, 63.9418923073, 44.5666425705, 44.5666425705, 27.272889691, 27.272889691},
 {58.1031793112, 58.1031793112, 43.1283065413, 43.1283065413, 32.7390076177, 32.7390076177},
 {53.2400113772, 53.2400113772, 40.8109126755, 40.8109126755, 30.2725048671, 30.2725048671},
 {49.1412215777, 49.1412215777, 38.7112014385, 38.7112014385, 29.6507376551, 29.6507376551},
 {45.6165488836, 45.6165488836, 36.6980118342, 36.6980118342, 28.6123054446, 28.6123054446},
 {42.5563710798, 42.5563710798, 34.8631723643, 34.8631723643, 27.7488594474, 27.7488594474},
 {39.8821236587, 39.8821236587, 33.2043584134, 33.2043584134, 26.9380420164, 26.9380420164},
 {37.5365654921, 37.5365654921, 31.7194345653, 31.7194345653, 26.2083695127, 26.2083695127},
 {35.4741284886, 35.4741284886, 30.3932289257, 30.3932289257, 25.5444626102, 25.5444626102}};

    TestHelper::dumpGolden("correctTemperatureError", results.temperature.errors);
    TestHelper::expectNear(correctTemperatureError, results.temperature.errors, 1e-6);
    TestHelper::dumpGolden("correctTemperatureSolution", results.temperature.values);
    TestHelper::expectNear(correctTemperatureSolution, results.temperature.values, 1e-6);

    // Checking number of iterations within subiterations
    EXPECT_EQ(progressMoisture.calls(), (std::vector<unsigned>{}));
    EXPECT_EQ(progressThermal.calls(), (std::vector<unsigned>{}));
}
