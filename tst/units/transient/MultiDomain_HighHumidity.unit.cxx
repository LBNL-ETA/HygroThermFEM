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

    const std::vector<double> correctHumidityError{5.17531008495e-11, 3.83820297076e-11, 2.9472658395e-11, 2.2334108018e-11, 1.68531563102e-11, 1.26500059674e-11, 9.45871872679e-12, 7.0587711966e-12, 5.25942912124e-12, 3.90372538567e-12};
    const std::vector<std::vector<double>> correctWaterContentSolution{{120.006901846, 120.006901846, 120.367262911, 120.367262911, 132.706612558, 132.706612558},
 {120.027041187, 120.027041187, 121.03321835, 121.03321835, 140.721285091, 140.721285091},
 {120.065162887, 120.065162887, 121.950679059, 121.950679059, 146.146712511, 146.146712511},
 {120.125348544, 120.125348544, 123.080509386, 123.080509386, 149.80865572, 149.80865572},
 {120.2101887, 120.2101887, 124.336142445, 124.336142445, 152.207125905, 152.207125905},
 {120.32101373, 120.32101373, 125.647111492, 125.647111492, 153.678646785, 153.678646785},
 {120.458374906, 120.458374906, 126.958206483, 126.958206483, 154.477859596, 154.477859596},
 {120.621887708, 120.621887708, 128.228716158, 128.228716158, 154.799595177, 154.799595177},
 {120.810393993, 120.810393993, 129.431522455, 129.431522455, 154.78806117, 154.788061169},
 {121.022243443, 121.022243443, 130.550684307, 130.550684307, 154.546838717, 154.546838717}};

    TestHelper::dumpGolden("correctHumidityError", results.moisture.errors);
    TestHelper::expectNear(correctHumidityError, results.moisture.errors, 1e-6);
    TestHelper::dumpGolden("correctWaterContentSolution", results.moisture.values);
    TestHelper::expectNear(correctWaterContentSolution, results.moisture.values, 1e-6);

    const std::vector<double> correctTemperatureError{6.44251932249e-07, 2.04929406944e-07, 1.30008732506e-07, 1.00725856001e-07, 8.22918117076e-08, 6.86307258569e-08, 5.78839119768e-08, 4.918518353e-08, 4.20267456945e-08, 3.60707476471e-08};
    const std::vector<std::vector<double>> correctTemperatureSolution{{0.799300504418, 0.799300504418, 2.60800534187, 2.60800534187, 7.95091536352, 7.95091536352},
 {2.10951781033, 2.10951781033, 5.0745564626, 5.0745564626, 10.7619945517, 10.7619945517},
 {3.63597938991, 3.63597938991, 7.09149829971, 7.09149829971, 12.2647216523, 12.2647216523},
 {5.20083005261, 5.20083005261, 8.74512309514, 8.74512309514, 13.3481312047, 13.3481312047},
 {6.71014579044, 6.71014579044, 10.1309797582, 10.1309797582, 14.219364247, 14.219364247},
 {8.11897365511, 8.11897365511, 11.3145789268, 11.3145789268, 14.9482972792, 14.9482972792},
 {9.40936199336, 9.40936199336, 12.3388534724, 12.3388534724, 15.5700284172, 15.5700284172},
 {10.5778138405, 10.5778138405, 13.232914936, 13.232914936, 16.1064565543, 16.1064565543},
 {11.6282939324, 11.6282939324, 14.0175876894, 14.0175876894, 16.5727161067, 16.5727161067},
 {12.5683666037, 12.5683666037, 14.7085988351, 14.7085988351, 16.9799566464, 16.9799566464}};

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

    const std::vector<double> correctHumidityError{1.24043141811e-11, 1.20989528296e-11, 9.85420655475e-12, 8.02101210212e-12, 6.68263540323e-12, 5.64742500241e-12, 4.81384337781e-12, 4.12822126464e-12, 3.55284451205e-12, 3.06250980035e-12};
    const std::vector<std::vector<double>> correctWaterContentSolution{{173.775711862, 173.775711862, 172.633351641, 172.633351641, 164.054650143, 164.054650143},
 {173.431879623, 173.431879623, 171.518251712, 171.518251712, 162.880271924, 162.880271924},
 {173.022082908, 173.022082908, 170.558513505, 170.558513505, 161.832767088, 161.832767088},
 {172.575084471, 172.575084471, 169.710984313, 169.710984313, 161.169387539, 161.169387539},
 {172.113140859, 172.113140859, 168.964082185, 168.964082185, 160.750376236, 160.750376236},
 {171.650583534, 171.650583534, 168.305136931, 168.305136931, 160.495443552, 160.495443552},
 {171.196407861, 171.196407861, 167.722730073, 167.722730073, 160.355431395, 160.355431395},
 {170.752976555, 170.752976555, 167.209023641, 167.209023641, 160.299083036, 160.299083036},
 {170.321964806, 170.321964806, 166.756683722, 166.756683722, 160.305162698, 160.305162698},
 {169.906441961, 169.906441961, 166.357948067, 166.357948067, 160.357984914, 160.357984914}};

    TestHelper::dumpGolden("correctHumidityError", results.moisture.errors);
    TestHelper::expectNear(correctHumidityError, results.moisture.errors, 1e-6);
    TestHelper::dumpGolden("correctWaterContentSolution", results.moisture.values);
    TestHelper::expectNear(correctWaterContentSolution, results.moisture.values, 1e-6);

    const std::vector<double> correctTemperatureError{1.894882554e-08, 2.39995097658e-08, 2.22409594205e-08, 2.05168319955e-08, 1.88487029839e-08, 1.72599726915e-08, 1.57648342838e-08, 1.4368453506e-08, 1.30727213829e-08, 1.18768446452e-08};
    const std::vector<std::vector<double>> correctTemperatureSolution{{29.5354697499, 29.5354697499, 28.2827561683, 28.2827561683, 24.1738393438, 24.1738393438},
 {28.8753624681, 28.8753624681, 27.1021603069, 27.1021603069, 23.8377635694, 23.8377635694},
 {28.1444608615, 28.1444608615, 26.1844164, 26.1844164, 23.3568471727, 23.3568471727},
 {27.4095713723, 27.4095713723, 25.4417197706, 25.4417197706, 22.9642884907, 22.9642884907},
 {26.7051336681, 26.7051336681, 24.8214230229, 24.8214230229, 22.6340216292, 22.6340216292},
 {26.047461279, 26.047461279, 24.2911714289, 24.2911714289, 22.3500934633, 22.3500934633},
 {25.4428551522, 25.4428551522, 23.83041293, 23.83041293, 22.1022337294, 22.1022337294},
 {24.892248635, 24.892248635, 23.4255825818, 23.4255825818, 21.8836064182, 21.8836064182},
 {24.393728116, 24.393728116, 23.0672614198, 23.0672614198, 21.6894366637, 21.6894366637},
 {23.9439755721, 23.9439755721, 22.748551477, 22.748551477, 21.5162107086, 21.5162107086}};

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

    const std::vector<double> correctHumidityError{1.68878475347e-08, 3.84346542168e-09, 2.94293811254e-10, 6.71711271483e-11, 1.34855041108e-10, 8.44265487656e-11, 5.37839649881e-11, 2.37605192365e-11, 2.27835908737e-10, 2.7784264776e-10};
    const std::vector<std::vector<double>> correctWaterContentSolution{{170.641312207, 170.641312207, 129.709323445, 129.709323445, 25.7909135045, 25.7909135045},
 {167.603299496, 167.603299496, 124.414593254, 124.414593254, 66.0882625371, 66.0882625371},
 {165.033239064, 165.033239064, 120.690657314, 120.690657314, 72.3473864704, 72.3473864704},
 {162.302554935, 162.302554935, 117.440353552, 117.440353552, 68.2352109425, 68.2352109425},
 {159.408338611, 159.408338611, 114.990212759, 114.990212759, 66.7277376209, 66.7277376209},
 {156.500672122, 156.500672122, 113.023310729, 113.023310729, 64.9832098512, 64.9832098512},
 {153.67173848, 153.67173848, 111.417607249, 111.417607249, 63.8053284571, 63.8053284571},
 {150.966132842, 150.966132842, 110.084662757, 110.084662757, 63.0033494636, 63.0033494636},
 {148.413559378, 148.413559378, 108.910993616, 108.910993616, 62.7159458317, 62.7159458317},
 {146.026221298, 146.026221298, 107.916759274, 107.916759274, 62.6769924082, 62.6769924082}};

    TestHelper::dumpGolden("correctHumidityError", results.moisture.errors);
    TestHelper::expectNear(correctHumidityError, results.moisture.errors, 1e-6);
    TestHelper::dumpGolden("correctWaterContentSolution", results.moisture.values);
    TestHelper::expectNear(correctWaterContentSolution, results.moisture.values, 1e-6);

    const std::vector<double> correctTemperatureError{9.53432428288e-07, 5.27389896785e-09, 8.42055948971e-08, 5.32274839728e-08, 5.58921121617e-08, 5.01270636157e-08, 4.68570821607e-08, 4.33004415417e-08, 4.00021999287e-08, 3.68358036195e-08};
    const std::vector<std::vector<double>> correctTemperatureSolution{{71.1495121697, 71.1495121697, 48.2705868338, 48.2705868338, -23.6784912828, -23.6784912828},
 {63.531186514, 63.5311865139, 44.4263968345, 44.4263968344, 27.5477931747, 27.5477931747},
 {57.6682184755, 57.6682184755, 42.9727036202, 42.9727036202, 32.660681303, 32.660681303},
 {52.8216807254, 52.8216807254, 40.6551626646, 40.6551626646, 30.2472339033, 30.2472339033},
 {48.7525610804, 48.7525610804, 38.5480123438, 38.5480123438, 29.5964530113, 29.5964530113},
 {45.2615631752, 45.2615631752, 36.5346061922, 36.5346061922, 28.5567628096, 28.5567628096},
 {42.2348738537, 42.2348738537, 34.7035339676, 34.7035339676, 27.6890834961, 27.6890834961},
 {39.592281574, 39.592281574, 33.0515957923, 33.0515957923, 26.8773046888, 26.8773046888},
 {37.2759944873, 37.2759944873, 31.5756070581, 31.5756070581, 26.1484789177, 26.1484789177},
 {35.2403387153, 35.2403387153, 30.2594927544, 30.2594927544, 25.4867393469, 25.4867393469}};

    TestHelper::dumpGolden("correctTemperatureError", results.temperature.errors);
    TestHelper::expectNear(correctTemperatureError, results.temperature.errors, 1e-6);
    TestHelper::dumpGolden("correctTemperatureSolution", results.temperature.values);
    TestHelper::expectNear(correctTemperatureSolution, results.temperature.values, 1e-6);

    // Checking number of iterations within subiterations
    EXPECT_EQ(progressMoisture.calls(), (std::vector<unsigned>{}));
    EXPECT_EQ(progressThermal.calls(), (std::vector<unsigned>{}));
}
