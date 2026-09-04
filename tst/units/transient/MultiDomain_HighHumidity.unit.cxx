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

    const std::vector<double> correctHumidityError{5.17529193826e-11, 3.83820297072e-11, 2.94726583953e-11, 2.23339265713e-11, 1.68533377533e-11, 1.26500059674e-11, 9.4587187267e-12, 7.05877119664e-12, 5.25924768659e-12, 3.90390681924e-12};
    const std::vector<std::vector<double>> correctWaterContentSolution{{120.006901847, 120.006901847, 120.367262988, 120.367262988, 132.706612324, 132.706612324},
 {120.02704118, 120.02704118, 121.033217446, 121.033217446, 140.721287818, 140.721287818},
 {120.065162886, 120.065162886, 121.950679983, 121.950679983, 146.146709744, 146.146709744},
 {120.125348552, 120.125348552, 123.080509657, 123.080509657, 149.808654892, 149.808654892},
 {120.210188691, 120.210188691, 124.336142749, 124.336142749, 152.207125013, 152.207125013},
 {120.321013724, 120.321013724, 125.647111065, 125.647111065, 153.678648084, 153.678648084},
 {120.458374889, 120.458374889, 126.9582057, 126.9582057, 154.477861982, 154.477861982},
 {120.621887658, 120.621887658, 128.228716728, 128.228716728, 154.799593569, 154.799593569},
 {120.810393938, 120.810393938, 129.431523371, 129.431523371, 154.788058535, 154.788058535},
 {121.022243422, 121.022243422, 130.55068373, 130.55068373, 154.546840495, 154.546840495}};

    TestHelper::dumpGolden("correctHumidityError", results.moisture.errors);
    TestHelper::expectNear(correctHumidityError, results.moisture.errors, 1e-6);
    TestHelper::dumpGolden("correctWaterContentSolution", results.moisture.values);
    TestHelper::expectNear(correctWaterContentSolution, results.moisture.values, 1e-6);

    const std::vector<double> correctTemperatureError{6.44251931968e-07, 2.04929407609e-07, 1.30008732339e-07, 1.00725855565e-07, 8.2291811857e-08, 6.86307261064e-08, 5.78839121981e-08, 4.91851835313e-08, 4.20267457021e-08, 3.60707476482e-08};
    const std::vector<std::vector<double>> correctTemperatureSolution{{0.799300504442, 0.799300504442, 2.6080053421, 2.6080053421, 7.95091536527, 7.95091536527},
 {2.10951781147, 2.10951781147, 5.0745564648, 5.0745564648, 10.7619945477, 10.7619945477},
 {3.63597938955, 3.63597938955, 7.09149829788, 7.09149829788, 12.2647216529, 12.2647216529},
 {5.2008300514, 5.2008300514, 8.74512309281, 8.74512309281, 13.3481312039, 13.3481312039},
 {6.7101457885, 6.7101457885, 10.1309797553, 10.1309797553, 14.2193642455, 14.2193642455},
 {8.11897365373, 8.11897365373, 11.314578926, 11.314578926, 14.948297279, 14.948297279},
 {9.40936199346, 9.40936199346, 12.3388534746, 12.3388534746, 15.5700284192, 15.5700284192},
 {10.5778138403, 10.5778138403, 13.2329149357, 13.2329149357, 16.1064565538, 16.1064565538},
 {11.628293931, 11.628293931, 14.0175876861, 14.0175876861, 16.5727161039, 16.5727161039},
 {12.5683666026, 12.5683666026, 14.7085988344, 14.7085988344, 16.9799566466, 16.9799566466}};

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

    const std::vector<double> correctHumidityError{1.24043141811e-11, 1.20987714987e-11, 9.85420655491e-12, 8.02119343774e-12, 6.68263540335e-12, 5.64742500254e-12, 4.81384337793e-12, 4.1282212647e-12, 3.55320719529e-12, 3.06250980039e-12};
    const std::vector<std::vector<double>> correctWaterContentSolution{{173.77571179, 173.77571179, 172.633351798, 172.633351798, 164.054649816, 164.054649816},
 {173.431881731, 173.431881731, 171.518251875, 171.518251875, 162.880267221, 162.880267221},
 {173.022085004, 173.022085004, 170.558513976, 170.558513976, 161.832761472, 161.832761472},
 {172.57508634, 172.57508634, 169.710984963, 169.710984964, 161.169381821, 161.169381821},
 {172.11314252, 172.11314252, 168.964082867, 168.964082867, 160.750370816, 160.750370816},
 {171.650584663, 171.650584663, 168.305138475, 168.305138475, 160.495436588, 160.495436588},
 {171.196408843, 171.196408843, 167.722731769, 167.722731769, 160.35542424, 160.35542424},
 {170.752978495, 170.752978495, 167.209024048, 167.209024048, 160.299077807, 160.299077807},
 {170.321966558, 170.321966558, 166.756684121, 166.756684121, 160.305157842, 160.305157842},
 {169.906443576, 169.906443576, 166.357948379, 166.357948379, 160.35798057, 160.35798057}};

    TestHelper::dumpGolden("correctHumidityError", results.moisture.errors);
    TestHelper::expectNear(correctHumidityError, results.moisture.errors, 1e-6);
    TestHelper::dumpGolden("correctWaterContentSolution", results.moisture.values);
    TestHelper::expectNear(correctWaterContentSolution, results.moisture.values, 1e-6);

    const std::vector<double> correctTemperatureError{1.89488253281e-08, 2.39995093219e-08, 2.22409591833e-08, 2.05168318638e-08, 1.8848702731e-08, 1.72599723126e-08, 1.57648342631e-08, 1.43684534869e-08, 1.30727213653e-08, 1.18768446294e-08};
    const std::vector<std::vector<double>> correctTemperatureSolution{{29.5354697501, 29.5354697501, 28.2827561695, 28.2827561695, 24.1738393439, 24.1738393439},
 {28.8753624756, 28.8753624756, 27.102160318, 27.102160318, 23.8377635817, 23.8377635817},
 {28.1444608767, 28.1444608767, 26.1844164185, 26.1844164185, 23.3568471875, 23.3568471875},
 {27.4095713939, 27.4095713939, 25.4417197946, 25.4417197946, 22.9642885085, 22.9642885085},
 {26.7051336949, 26.7051336949, 24.8214230505, 24.8214230505, 22.6340216484, 22.6340216484},
 {26.0474613093, 26.0474613093, 24.2911714611, 24.2911714611, 22.3500934868, 22.3500934868},
 {25.4428551853, 25.4428551853, 23.8304129654, 23.8304129654, 22.1022337545, 22.1022337545},
 {24.8922486715, 24.8922486715, 23.4255826162, 23.4255826162, 21.8836064399, 21.8836064399},
 {24.393728154, 24.393728154, 23.0672614531, 23.0672614531, 21.6894366844, 21.6894366844},
 {23.9439756106, 23.9439756106, 22.748551509, 22.748551509, 21.5162107281, 21.5162107281}};

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

    const std::vector<double> correctHumidityError{1.68873659622e-08, 3.843260157e-09, 2.9428218136e-10, 6.71090942008e-11, 1.34853588528e-10, 8.44261866195e-11, 5.37837841075e-11, 2.37605197429e-11, 2.28025656291e-10, 2.77854314818e-10};
    const std::vector<std::vector<double>> correctWaterContentSolution{{170.641308484, 170.641308484, 129.709277291, 129.709277291, 25.790836353, 25.790836353},
 {167.603298043, 167.603298043, 124.41458293, 124.41458293, 66.0878591567, 66.0878591567},
 {165.033240299, 165.033240299, 120.690600942, 120.690600942, 72.3470762537, 72.3470762537},
 {162.302543377, 162.302543377, 117.440290779, 117.440290779, 68.2349387249, 68.2349387249},
 {159.408316063, 159.408316063, 114.990148134, 114.990148134, 66.7274898156, 66.7274898156},
 {156.500640526, 156.500640526, 113.023246976, 113.023246976, 64.9829746586, 64.9829746586},
 {153.671699976, 153.671699976, 111.417545372, 111.417545372, 63.8050999543, 63.8050999543},
 {150.966090833, 150.966090833, 110.084602397, 110.084602397, 63.0031224316, 63.0031224316},
 {148.413512819, 148.413512819, 108.910909778, 108.910909778, 62.7158004532, 62.7158004532},
 {146.02617327, 146.02617327, 107.916667789, 107.916667789, 62.6768719438, 62.6768719438}};

    TestHelper::dumpGolden("correctHumidityError", results.moisture.errors);
    TestHelper::expectNear(correctHumidityError, results.moisture.errors, 1e-6);
    TestHelper::dumpGolden("correctWaterContentSolution", results.moisture.values);
    TestHelper::expectNear(correctWaterContentSolution, results.moisture.values, 1e-6);

    const std::vector<double> correctTemperatureError{9.53432762926e-07, 5.27392366883e-09, 8.42055712476e-08, 5.32274520461e-08, 5.58921071562e-08, 5.01270571023e-08, 4.68570808412e-08, 4.33004410636e-08, 4.00021964937e-08, 3.68357979512e-08};
    const std::vector<std::vector<double>> correctTemperatureSolution{{71.1495113653, 71.1495113653, 48.2705844074, 48.2705844074, -23.6784934921, -23.6784934921},
 {63.5311883302, 63.5311883302, 44.4264055175, 44.4264055175, 27.5478191465, 27.5478191465},
 {57.6682218188, 57.6682218188, 42.9727108856, 42.9727108856, 32.660685147, 32.660685147},
 {52.8216848394, 52.8216848394, 40.655168294, 40.655168294, 30.2472367713, 30.2472367713},
 {48.7525653139, 48.7525653139, 38.5480168807, 38.5480168807, 29.5964557218, 29.5964557218},
 {45.2615671135, 45.2615671135, 36.5346097492, 36.5346097492, 28.5567647543, 28.5567647543},
 {42.2348773118, 42.2348773118, 34.7035367914, 34.7035367914, 27.6890851451, 27.6890851451},
 {39.5922845065, 39.5922845065, 33.0515980337, 33.0515980337, 26.8773060325, 26.8773060325},
 {37.2759972064, 37.2759972064, 31.5756097628, 31.5756097628, 26.1484818582, 26.1484818582},
 {35.2403412286, 35.2403412286, 30.2594951996, 30.2594951996, 25.4867414106, 25.4867414106}};

    TestHelper::dumpGolden("correctTemperatureError", results.temperature.errors);
    TestHelper::expectNear(correctTemperatureError, results.temperature.errors, 1e-6);
    TestHelper::dumpGolden("correctTemperatureSolution", results.temperature.values);
    TestHelper::expectNear(correctTemperatureSolution, results.temperature.values, 1e-6);

    // Checking number of iterations within subiterations
    EXPECT_EQ(progressMoisture.calls(), (std::vector<unsigned>{}));
    EXPECT_EQ(progressThermal.calls(), (std::vector<unsigned>{}));
}
