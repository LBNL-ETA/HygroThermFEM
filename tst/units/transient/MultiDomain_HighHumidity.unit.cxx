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

    const std::vector<double> correctHumidityError{1.68875975813e-08, 3.84351734048e-09, 2.94285814476e-10, 6.71658520575e-11, 1.34853040462e-10, 8.44256392478e-11, 5.37830550996e-11, 2.37603372955e-11, 2.27889081203e-10, 2.77846474497e-10};
    const std::vector<std::vector<double>> correctWaterContentSolution{{170.641312207, 170.641312207, 129.709323445, 129.709323445, 25.7909135045, 25.7909135045},
 {167.603306453, 167.603306453, 124.414663466, 124.414663466, 66.0880381124, 66.0880381124},
 {165.033248676, 165.033248676, 120.690699423, 120.690699423, 72.3472172437, 72.3472172437},
 {162.302568476, 162.302568476, 117.440372504, 117.440372504, 68.235097766, 68.235097766},
 {159.408353366, 159.408353366, 114.990220855, 114.990220855, 66.7276499848, 66.7276499848},
 {156.500686451, 156.500686451, 113.023313095, 113.023313095, 64.983136371, 64.983136371},
 {153.671751744, 153.671751744, 111.417606474, 111.417606474, 63.8052635444, 63.8052635444},
 {150.966144852, 150.966144852, 110.084660128, 110.084660128, 63.0032901034, 63.0032901034},
 {148.413569365, 148.413569365, 108.910983952, 108.910983952, 62.7159101911, 62.7159101911},
 {146.026229357, 146.026229357, 107.916747047, 107.916747047, 62.6769662326, 62.6769662325}};

    TestHelper::dumpGolden("correctHumidityError", results.moisture.errors);
    TestHelper::expectNear(correctHumidityError, results.moisture.errors, 1e-6);
    TestHelper::dumpGolden("correctWaterContentSolution", results.moisture.values);
    TestHelper::expectNear(correctWaterContentSolution, results.moisture.values, 1e-6);

    const std::vector<double> correctTemperatureError{9.53432428288e-07, 5.27390195932e-09, 8.42055743524e-08, 5.32274644429e-08, 5.58921052873e-08, 5.01270585754e-08, 4.68570806599e-08, 4.33004413767e-08, 4.00021999384e-08, 3.68358039103e-08};
    const std::vector<std::vector<double>> correctTemperatureSolution{{71.1495121697, 71.1495121697, 48.2705868338, 48.2705868338, -23.6784912828, -23.6784912828},
 {63.531188631, 63.531188631, 44.4264036094, 44.4264036094, 27.5478087059, 27.5478087059},
 {57.6682218848, 57.6682218848, 42.9727096129, 42.9727096129, 32.6606845678, 32.6606845678},
 {52.8216847623, 52.8216847623, 40.6551679071, 40.6551679071, 30.2472368307, 30.2472368307},
 {48.7525653803, 48.7525653803, 38.5480168693, 38.5480168693, 29.5964555451, 29.5964555451},
 {45.2615674826, 45.2615674826, 36.5346100835, 36.5346100835, 28.5567648792, 28.5567648792},
 {42.2348780074, 42.2348780074, 34.70353736, 34.70353736, 27.6890853201, 27.6890853201},
 {39.5922854808, 39.5922854808, 33.0515987749, 33.0515987749, 26.8773062908, 26.8773062908},
 {37.2759981792, 37.2759981792, 31.5756099338, 31.5756099338, 26.1484807999, 26.1484807999},
 {35.2403421526, 35.2403421526, 30.2594953419, 30.2594953419, 25.4867408603, 25.4867408603}};

    TestHelper::dumpGolden("correctTemperatureError", results.temperature.errors);
    TestHelper::expectNear(correctTemperatureError, results.temperature.errors, 1e-6);
    TestHelper::dumpGolden("correctTemperatureSolution", results.temperature.values);
    TestHelper::expectNear(correctTemperatureSolution, results.temperature.values, 1e-6);

    // Checking number of iterations within subiterations
    EXPECT_EQ(progressMoisture.calls(), (std::vector<unsigned>{}));
    EXPECT_EQ(progressThermal.calls(), (std::vector<unsigned>{}));
}
