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

    const std::vector<double> correctHumidityError{2.258469e-11, 2.322582e-11, 2.231239e-11, 2.262864e-11, 1.645984e-11, 1.733963e-11, 1.295150e-11, 1.378937e-11, 1.025800e-11, 1.088118e-11};
    const std::vector<std::vector<double>> correctWaterContentSolution{
      {120.214619, 120.214619, 120.394521, 120.394521, 120.608799, 120.608799},
      {120.338678, 120.338678, 120.486319, 120.486319, 120.641100, 120.641100},
      {120.503448, 120.503448, 120.601687, 120.601687, 120.722768, 120.722768},
      {120.703536, 120.703536, 120.789481, 120.789481, 120.883126, 120.883126},
      {120.724042, 120.724042, 120.806993, 120.806993, 120.906566, 120.906566},
      {120.813198, 120.813198, 120.884634, 120.884634, 120.961029, 120.961029},
      {120.830219, 120.830219, 120.898469, 120.898469, 120.977725, 120.977725},
      {120.903197, 120.903197, 120.958845, 120.958845, 121.018848, 121.018848},
      {120.916487, 120.916487, 120.969821, 120.969821, 121.032337, 121.032337},
      {120.973684, 120.973684, 121.017799, 121.017799, 121.065295, 121.065295}};

    TestHelper::expectNear(correctHumidityError, results.moisture.errors, 1e-6);
    TestHelper::expectNear(correctWaterContentSolution, results.moisture.values, 1e-6);

    const std::vector<double> correctTemperatureError{6.314878e-07, 1.970516e-07, 1.253964e-07, 9.751431e-08, 7.984890e-08, 6.668574e-08, 5.631733e-08, 4.789961e-08, 4.096431e-08, 3.517708e-08};
    const std::vector<std::vector<double>> correctTemperatureSolution{
      {0.809876, 0.809876, 2.644933, 2.644933, 8.059481, 8.059481},
      {2.133188, 2.133188, 5.132997, 5.132997, 10.864948, 10.864948},
      {3.671155, 3.671155, 7.159582, 7.159582, 12.352463, 12.352463},
      {5.245095, 5.245095, 8.817682, 8.817682, 13.427511, 13.427511},
      {6.762352, 6.762352, 10.206564, 10.206564, 14.293930, 14.293930},
      {8.177890, 8.177890, 11.392225, 11.392225, 15.019141, 15.019141},
      {9.474261, 9.474261, 12.418194, 12.418194, 15.637973, 15.637973},
      {10.647815, 10.647815, 13.313537, 13.313537, 16.172017, 16.172017},
      {11.702779, 11.702779, 14.099249, 14.099249, 16.636355, 16.636355},
      {12.646561, 12.646561, 14.790907, 14.790907, 17.041890, 17.041890}};

    TestHelper::expectNear(correctTemperatureError, results.temperature.errors, 1e-6);
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

    const std::vector<double> correctHumidityError{6.071259e-10, 2.366881e-11, 3.948503e-11, 4.818683e-11, 3.936627e-11, 2.912424e-11, 2.692081e-11, 3.116192e-11, 2.253094e-11, 1.884293e-11};
    const std::vector<std::vector<double>> correctWaterContentSolution{
      {101.650289, 101.650289, 96.087275, 96.087275, 91.936622, 91.936622},
      {97.840964, 97.840964, 96.124041, 96.124041, 95.262493, 95.262493},
      {96.492198, 96.492198, 95.823509, 95.823509, 95.289098, 95.289098},
      {95.870839, 95.870839, 95.485076, 95.485076, 95.117559, 95.117559},
      {95.415139, 95.415139, 95.149303, 95.149303, 94.830863, 94.830863},
      {95.041234, 95.041234, 94.825110, 94.825110, 94.528052, 94.528052},
      {94.776780, 94.776780, 94.589143, 94.589143, 94.339786, 94.339786},
      {94.570722, 94.570722, 94.396682, 94.396682, 94.216143, 94.216143},
      {94.244013, 94.244013, 94.111271, 94.111271, 93.935107, 93.935107},
      {94.133276, 94.133276, 94.012349, 94.012349, 93.855730, 93.855730}};

    TestHelper::expectNear(correctHumidityError, results.moisture.errors, 1e-6);
    TestHelper::expectNear(correctWaterContentSolution, results.moisture.values, 1e-6);

    const std::vector<double> correctTemperatureError{8.107215e-09, 2.935156e-08, 2.516332e-08, 2.272577e-08, 2.028107e-08, 1.801523e-08, 1.594434e-08, 1.407456e-08, 1.240238e-08, 1.090569e-08};
    const std::vector<std::vector<double>> correctTemperatureSolution{
      {29.266227, 29.266227, 27.717828, 27.717828, 23.388971, 23.388971},
      {28.349757, 28.349757, 26.438542, 26.438542, 23.491167, 23.491167},
      {27.407023, 27.407023, 25.449689, 25.449689, 22.934090, 22.934090},
      {26.515300, 26.515300, 24.667888, 24.667888, 22.530777, 22.530777},
      {25.704188, 25.704188, 24.026567, 24.026567, 22.193728, 22.193728},
      {24.981352, 24.981352, 23.488405, 23.488405, 21.908855, 21.908855},
      {24.344467, 24.344467, 23.030354, 23.030354, 21.664882, 21.664882},
      {23.786803, 23.786803, 22.637053, 22.637053, 21.454350, 21.454350},
      {23.299951, 23.299951, 22.297429, 22.297429, 21.271879, 21.271879},
      {22.875999, 22.875999, 22.003370, 22.003370, 21.113282, 21.113282}};

    TestHelper::expectNear(correctTemperatureError, results.temperature.errors, 1e-6);
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

    const std::vector<double> correctHumidityError{4.571666e-08, 1.495656e-08, 3.977508e-09, 1.149587e-09, 2.043850e-10, 3.180838e-11, 1.121104e-10, 1.290410e-10, 1.255300e-10, 1.148279e-10};
    const std::vector<std::vector<double>> correctWaterContentSolution{
      {94.067191, 94.067191, 78.841515, 78.841515, 15.003698, 15.003698},
      {77.171092, 77.170677, 63.609643, 63.609443, 27.228134, 27.228121},
      {70.670721, 70.670630, 60.214856, 60.214843, 46.761334, 46.761369},
      {65.165722, 65.165638, 58.724426, 58.724413, 51.215218, 51.215223},
      {62.385014, 62.385007, 57.944295, 57.944291, 53.224124, 53.224126},
      {61.375145, 61.375142, 57.476338, 57.476337, 54.037252, 54.037253},
      {60.494227, 60.494226, 57.154797, 57.154797, 54.470938, 54.470938},
      {59.735444, 59.735443, 56.906731, 56.906731, 54.729504, 54.729504},
      {59.086618, 59.086618, 56.703990, 56.703990, 54.906435, 54.906435},
      {58.534404, 58.534404, 56.533549, 56.533549, 55.036971, 55.036971}};

    TestHelper::expectNear(correctHumidityError, results.moisture.errors, 1e-6);
    TestHelper::expectNear(correctWaterContentSolution, results.moisture.values, 1e-6);

    const std::vector<double> correctTemperatureError{1.210316e-06, 2.305046e-08, 9.875415e-08, 6.578440e-08, 6.382298e-08, 5.538421e-08, 4.932550e-08, 4.357260e-08, 3.840565e-08, 3.372338e-08};
    const std::vector<std::vector<double>> correctTemperatureSolution{
      {67.554263, 67.554263, 42.195367, 42.195367, -25.747790, -25.747790},
      {57.873437, 57.873435, 39.366726, 39.366726, 26.197255, 26.197256},
      {51.000498, 51.000497, 38.183576, 38.183576, 30.166558, 30.166557},
      {45.666114, 45.666113, 35.920427, 35.920427, 28.011680, 28.011680},
      {41.452159, 41.452158, 33.834519, 33.834519, 27.239000, 27.239000},
      {38.040198, 38.040198, 31.897744, 31.897744, 26.244943, 26.244943},
      {35.232856, 35.232856, 30.196674, 30.196674, 25.414864, 25.414864},
      {32.898411, 32.898411, 28.723167, 28.723167, 24.677240, 24.677240},
      {30.944546, 30.944546, 27.458938, 27.458938, 24.037794, 24.037794},
      {29.302317, 29.302317, 26.379025, 26.379025, 23.485058, 23.485058}};

    TestHelper::expectNear(correctTemperatureError, results.temperature.errors, 1e-6);
    TestHelper::expectNear(correctTemperatureSolution, results.temperature.values, 1e-6);

    // Checking number of iterations within subiterations
    EXPECT_EQ(progressMoisture.calls(), (std::vector<unsigned>{}));
    EXPECT_EQ(progressThermal.calls(), (std::vector<unsigned>{}));
}
