#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "SlabCreator.hxx"
#include "TestHelpers.hxx"
#include "TestMaterials.hxx"

using HygroThermFEM::Nodes;
using HygroThermFEM::SimulationProperties;

class MultiDomain_2D_ThermalConductivityDependent_1 : public testing::Test
{
protected:
    void SetUp() override
    {
        constexpr auto relaxationParameter{0.8};
        constexpr auto errorTolerance{1e-5};
        constexpr auto numberOfIterations{20u};
        SimulationProperties::Instance().setIterationParameters(
          relaxationParameter, errorTolerance, numberOfIterations);
    }

    void TearDown() override
    {
        SimulationProperties::Instance().reset();
    }
};

TEST_F(MultiDomain_2D_ThermalConductivityDependent_1, TestExample_1)
{
    constexpr auto excludeWaterLiquidTransportation{false};
    constexpr auto excludeHeatOfEvaporation{false};
    constexpr auto excludeCapillaryConduction{false};
    constexpr auto excludeVaporDiffusionConduction{false};
    constexpr auto thermalConductivityMoistureAndTemperatureDependent{true};

    SimulationProperties::Instance().setCalculationParameters(
      excludeWaterLiquidTransportation,
      excludeHeatOfEvaporation,
      excludeCapillaryConduction,
      excludeVaporDiffusionConduction,
      thermalConductivityMoistureAndTemperatureDependent);

    HygroThermFEM::MultiDomain multiDomain;

    auto params = TestHelper::CottaerSandstone();
    params.thermalConductivityMoistureDependent = {{0.0, 1.8}, {180, 2.5}};
    params.thermalConductivityTemperatureDependent = {{0.0, 1.8}, {1, 3.1}};
    const auto & material = multiDomain.materials().createSolidMaterial(params);

    constexpr HygroThermFEM::State state({
        .temperature = 0.0,
        .humidity = 0.99,
        .pressure = 101325.0,
        .liquidPercent = 1.0
    });

    TestHelper::SlabBuilder(multiDomain)
        .gridXCoordinates({0, 0.05, 0.1})
        .height(0.05)
        .material(material.name())
        .state(state)
        .startCorner(TestHelper::StartCorner::BottomRight)
        .direction(TestHelper::Direction::CounterClockwise)
        .build();

    /// Create Boundary Conditions
    constexpr auto hc = 5.0;
    constexpr auto airTemperature = 10.0;
    constexpr auto humidity = 0.0;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff{airTemperature, hc, humidity};

    multiDomain.createBC_FixedHc(1, 2, bcCoeff);
    multiDomain.createBC_FixedHc(5, 6, bcCoeff);

    constexpr auto dTime = 3600;
    constexpr auto nSteps = 24;

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

    const std::vector<std::vector<double>> correctWaterContentSolution{
      {46.885036, 46.885036, 47.305674, 47.305674, 46.885036, 46.885036},
      {50.397795, 50.397795, 48.911762, 48.911762, 50.397795, 50.397795},
      {48.755896, 48.755896, 49.080708, 49.080708, 48.755896, 48.755896},
      {47.866724, 47.866724, 48.357078, 48.357078, 47.866724, 47.866724},
      {46.970685, 46.970685, 47.495224, 47.495224, 46.970685, 46.970685},
      {45.983809, 45.983809, 46.523640, 46.523640, 45.983809, 45.983809},
      {44.676292, 44.676292, 45.167777, 45.167777, 44.676292, 44.676292},
      {49.023010, 49.023010, 48.758598, 48.758598, 49.023010, 49.023010},
      {50.631184, 50.631184, 50.798793, 50.798793, 50.631184, 50.631184},
      {49.783028, 49.783028, 50.241127, 50.241127, 49.783028, 49.783028},
      {49.019897, 49.019897, 49.517548, 49.517548, 49.019897, 49.019897},
      {48.207729, 48.207729, 48.725098, 48.725098, 48.207729, 48.207729},
      {47.350196, 47.350196, 47.867247, 47.867247, 47.350196, 47.350196},
      {48.624259, 48.624259, 48.898061, 48.898061, 48.624259, 48.624259},
      {47.680008, 47.680008, 48.175189, 48.175189, 47.680008, 47.680008},
      {46.753273, 46.753273, 47.288168, 47.288168, 46.753273, 46.753273},
      {45.745723, 45.745723, 46.286810, 46.286810, 45.745723, 45.745723},
      {44.598831, 44.598831, 45.175582, 45.175582, 44.598831, 44.598831},
      {43.034735, 43.034735, 44.123832, 44.123832, 43.034735, 43.034735},
      {41.081889, 41.081889, 43.248167, 43.248167, 41.081889, 41.081889},
      {39.312049, 39.312049, 42.106711, 42.106711, 39.312049, 39.312049},
      {37.620245, 37.620245, 40.797143, 40.797143, 37.620245, 37.620245},
      {35.947220, 35.947220, 39.366294, 39.366294, 35.947220, 35.947220},
      {34.254060, 34.254060, 37.839858, 37.839858, 34.254060, 34.254060}};

    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    const std::vector<std::vector<double>> correctTemperatureSolution{
      {0.142456, 0.142456, 0.097914, 0.097914, 0.142456, 0.142456},
      {0.226626, 0.226626, 0.194314, 0.194314, 0.226626, 0.226626},
      {0.293843, 0.293843, 0.270713, 0.270713, 0.293843, 0.293843},
      {0.348760, 0.348760, 0.331197, 0.331197, 0.348760, 0.348760},
      {0.393059, 0.393059, 0.379505, 0.379505, 0.393059, 0.393059},
      {0.429097, 0.429097, 0.418369, 0.418369, 0.429097, 0.429097},
      {0.458901, 0.458901, 0.450513, 0.450513, 0.458901, 0.458901},
      {0.476908, 0.476908, 0.471967, 0.471967, 0.476908, 0.476908},
      {0.489366, 0.489366, 0.486276, 0.486276, 0.489366, 0.489366},
      {0.498301, 0.498301, 0.495760, 0.495760, 0.498301, 0.498301},
      {0.506439, 0.506439, 0.504168, 0.504168, 0.506439, 0.506439},
      {0.514006, 0.514006, 0.511900, 0.511900, 0.514006, 0.514006},
      {0.520947, 0.520947, 0.519025, 0.519025, 0.520947, 0.520947},
      {0.525880, 0.525880, 0.524833, 0.524833, 0.525880, 0.525880},
      {0.529610, 0.529610, 0.528543, 0.528543, 0.529610, 0.529610},
      {0.533886, 0.533886, 0.532700, 0.532700, 0.533886, 0.533886},
      {0.538606, 0.538606, 0.537305, 0.537305, 0.538606, 0.538606},
      {0.543755, 0.543755, 0.542336, 0.542336, 0.543755, 0.543755},
      {0.550510, 0.550510, 0.548747, 0.548747, 0.550510, 0.550510},
      {0.558750, 0.558750, 0.556639, 0.556639, 0.558750, 0.558750},
      {0.567627, 0.567627, 0.565348, 0.565348, 0.567627, 0.567627},
      {0.576931, 0.576931, 0.574559, 0.574559, 0.576931, 0.576931},
      {0.586559, 0.586559, 0.584129, 0.584129, 0.586559, 0.586559},
      {0.596481, 0.596481, 0.594009, 0.594009, 0.596481, 0.596481}};

    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-6);
}
