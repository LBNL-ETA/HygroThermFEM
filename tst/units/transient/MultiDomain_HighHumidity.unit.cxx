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

    const std::vector<double> correctHumidityError{5.17478383e-11, 3.83635212e-11, 2.94503405e-11, 2.23139683e-11, 1.68364645e-11, 1.26373061e-11, 9.4494662e-12, 7.05260304e-12, 5.25561952e-12, 3.90227433e-12};
    const std::vector<std::vector<double>> correctWaterContentSolution{{120.006902, 120.006902, 120.367263, 120.367263, 132.706613, 132.706613},
 {120.027041, 120.027041, 121.033183, 121.033183, 140.720431, 140.720431},
 {120.065157, 120.065157, 121.950424, 121.950424, 146.143133, 146.143133},
 {120.125329, 120.125329, 123.079778, 123.079778, 149.802087, 149.802087},
 {120.210141, 120.210141, 124.334743, 124.334743, 152.198027, 152.198027},
 {120.320922, 120.320922, 125.644892, 125.644892, 153.66773, 153.66773},
 {120.458217, 120.458217, 126.955091, 126.955091, 154.465777, 154.465777},
 {120.621643, 120.621643, 128.224699, 128.224699, 154.786869, 154.786869},
 {120.810041, 120.810041, 129.426653, 129.426653, 154.775102, 154.775102},
 {121.021762, 121.021762, 130.545043, 130.545043, 154.533968, 154.533968}};

    TestHelper::dumpGolden("correctHumidityError", results.moisture.errors);
    TestHelper::expectNear(correctHumidityError, results.moisture.errors, 1e-6);
    TestHelper::dumpGolden("correctWaterContentSolution", results.moisture.values);
    TestHelper::expectNear(correctWaterContentSolution, results.moisture.values, 1e-6);

    const std::vector<double> correctTemperatureError{6.4522847e-07, 2.0509422e-07, 1.29942904e-07, 1.00626399e-07, 8.22002844e-08, 6.85536417e-08, 5.78210798e-08, 4.91351013e-08, 4.19877101e-08, 3.60411294e-08};
    const std::vector<std::vector<double>> correctTemperatureSolution{{0.798775906, 0.798775906, 2.60752793, 2.60752793, 7.952544, 7.952544},
 {2.10776735, 2.10776735, 5.07323751, 5.07323751, 10.7669095, 10.7669095},
 {3.63269828, 3.63269828, 7.08959472, 7.08959472, 12.2706537, 12.2706537},
 {5.19599656, 5.19599656, 8.74276519, 8.74276519, 13.3539079, 13.3539079},
 {6.70388745, 6.70388745, 10.128199, 10.128199, 14.2246145, 14.2246145},
 {8.11148648, 8.11148648, 11.3113815, 11.3113815, 14.9529079, 14.9529079},
 {9.40086612, 9.40086612, 12.3352517, 12.3352517, 15.5739686, 15.5739686},
 {10.5685296, 10.5685296, 13.2289353, 13.2289353, 16.1097318, 16.1097318},
 {11.6184295, 11.6184295, 14.0132708, 14.0132708, 16.5753539, 16.5753539},
 {12.5581117, 12.5581117, 14.7039961, 14.7039961, 16.9819995, 16.9819995}};

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

    const std::vector<double> correctHumidityError{1.23513663e-11, 1.20247879e-11, 9.81775745e-12, 8.00160806e-12, 6.67374879e-12, 5.64579197e-12, 4.81692529e-12, 4.13474879e-12, 3.56191097e-12, 3.07266445e-12};
    const std::vector<std::vector<double>> correctWaterContentSolution{{173.775712, 173.775712, 172.633352, 172.633352, 164.05465, 164.05465},
 {173.431963, 173.431963, 171.518891, 171.518891, 162.887685, 162.887685},
 {173.02242, 173.02242, 170.560646, 170.560646, 161.849852, 161.849852},
 {172.575796, 172.575796, 169.714674, 169.714674, 161.190862, 161.190862},
 {172.114295, 172.114295, 168.96918, 168.96918, 160.773699, 160.773699},
 {171.652209, 171.652209, 168.311424, 168.311424, 160.519092, 160.519092},
 {171.198502, 171.198502, 167.729971, 167.729971, 160.378442, 160.378442},
 {170.755556, 170.755556, 167.21697, 167.21697, 160.320827, 160.320827},
 {170.325001, 170.325001, 166.765125, 166.765125, 160.325248, 160.325248},
 {169.909898, 169.909898, 166.366695, 166.366695, 160.37619, 160.37619}};

    TestHelper::dumpGolden("correctHumidityError", results.moisture.errors);
    TestHelper::expectNear(correctHumidityError, results.moisture.errors, 1e-6);
    TestHelper::dumpGolden("correctWaterContentSolution", results.moisture.values);
    TestHelper::expectNear(correctWaterContentSolution, results.moisture.values, 1e-6);

    const std::vector<double> correctTemperatureError{1.9019259e-08, 2.3867777e-08, 2.2147929e-08, 2.04453472e-08, 1.87949127e-08, 1.7220615e-08, 1.57371876e-08, 1.43502965e-08, 1.30622113e-08, 1.18724793e-08};
    const std::vector<std::vector<double>> correctTemperatureSolution{{29.5372986, 29.5372986, 28.2847356, 28.2847356, 24.1669619, 24.1669619},
 {28.881341, 28.881341, 27.1069531, 27.1069531, 23.8269724, 23.8269724},
 {28.1544133, 28.1544133, 26.1903473, 26.1903473, 23.3494517, 23.3494517},
 {27.422632, 27.422632, 25.4483064, 25.4483064, 22.9589365, 22.9589365},
 {26.7203523, 26.7203523, 24.8285051, 24.8285051, 22.630241, 22.630241},
 {26.0640162, 26.0640162, 24.2986699, 24.2986699, 22.3475832, 22.3475832},
 {25.4601022, 25.4601022, 23.8382515, 23.8382515, 22.1007692, 22.1007692},
 {24.9097077, 24.9097077, 23.4336713, 23.4336713, 21.8830022, 21.8830022},
 {24.4110538, 24.4110538, 23.075502, 23.075502, 21.6895332, 21.6895332},
 {23.9609251, 23.9609251, 22.7568457, 22.7568457, 21.5168689, 21.5168689}};

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

    const std::vector<double> correctHumidityError{1.68964184e-08, 4.03001658e-09, 3.18288146e-10, 6.97290545e-11, 1.42441618e-10, 9.16290735e-11, 6.08170742e-11, 3.03779118e-11, 1.58943833e-10, 2.62557039e-10};
    const std::vector<std::vector<double>> correctWaterContentSolution{{170.641312, 170.641312, 129.709323, 129.709323, 25.7909135, 25.7909135},
 {167.602123, 167.602123, 124.423515, 124.423515, 66.1412532, 66.1412532},
 {165.031891, 165.031891, 120.77303, 120.77303, 72.7810058, 72.7810058},
 {162.31747, 162.31747, 117.528851, 117.528851, 68.5485762, 68.5485762},
 {159.436005, 159.436005, 115.077195, 115.077195, 67.0062269, 67.0062269},
 {156.536687, 156.536687, 113.1027, 113.1027, 65.1981742, 65.1981742},
 {153.712471, 153.712471, 111.486312, 111.486312, 63.9537764, 63.9537764},
 {151.00881, 151.00881, 110.141125, 110.141125, 63.0826026, 63.0826026},
 {148.455829, 148.455829, 108.956082, 108.956082, 62.7235657, 62.7235657},
 {146.066174, 146.066174, 107.943175, 107.943175, 62.640355, 62.640355}};

    TestHelper::dumpGolden("correctHumidityError", results.moisture.errors);
    TestHelper::expectNear(correctHumidityError, results.moisture.errors, 1e-6);
    TestHelper::dumpGolden("correctWaterContentSolution", results.moisture.values);
    TestHelper::expectNear(correctWaterContentSolution, results.moisture.values, 1e-6);

    const std::vector<double> correctTemperatureError{9.822004e-07, 1.66435058e-09, 8.41482712e-08, 5.29825674e-08, 5.59020527e-08, 5.02374814e-08, 4.70444797e-08, 4.3543668e-08, 4.02855958e-08, 3.71463727e-08};
    const std::vector<std::vector<double>> correctTemperatureSolution{{71.4287883, 71.4287883, 48.5876104, 48.5876104, -24.7690146, -24.7690146},
 {64.0327887, 64.0327887, 44.5541365, 44.5541365, 27.1700736, 27.1700736},
 {58.2364873, 58.2364873, 43.1221339, 43.1221339, 32.6931379, 32.6931379},
 {53.3835522, 53.3835522, 40.8195578, 40.8195578, 30.2535475, 30.2535475},
 {49.2825638, 49.2825638, 38.7340475, 38.7340475, 29.6346956, 29.6346956},
 {45.7506545, 45.7506545, 36.7323254, 36.7323254, 28.608828, 28.608828},
 {42.6814089, 42.6814089, 34.9046894, 34.9046894, 27.7516189, 27.7516189},
 {39.9975241, 39.9975241, 33.2495017, 33.2495017, 26.9454445, 26.9454445},
 {37.6423805, 37.6423805, 31.7657207, 31.7657207, 26.2187792, 26.2187792},
 {35.570798, 35.570798, 30.439293, 30.439293, 25.5573648, 25.5573648}};

    TestHelper::dumpGolden("correctTemperatureError", results.temperature.errors);
    TestHelper::expectNear(correctTemperatureError, results.temperature.errors, 1e-6);
    TestHelper::dumpGolden("correctTemperatureSolution", results.temperature.values);
    TestHelper::expectNear(correctTemperatureSolution, results.temperature.values, 1e-6);

    // Checking number of iterations within subiterations
    EXPECT_EQ(progressMoisture.calls(), (std::vector<unsigned>{}));
    EXPECT_EQ(progressThermal.calls(), (std::vector<unsigned>{}));
}
