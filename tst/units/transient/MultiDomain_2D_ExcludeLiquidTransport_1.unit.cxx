#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "SlabCreator.hxx"
#include "TestHelpers.hxx"
#include "TestMaterials.hxx"

using HygroThermFEM::Nodes;
using HygroThermFEM::SimulationProperties;

class MultiDomain_2D_ExcludeLiquidTransport_1 : public testing::Test
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

TEST_F(MultiDomain_2D_ExcludeLiquidTransport_1, TestExample_1)
{
    constexpr auto excludeWaterLiquidTransportation{true};
    constexpr auto excludeHeatOfEvaporation{false};
    constexpr auto excludeCapillaryConduction{false};
    constexpr auto excludeVaporDiffusionConduction{false};
    constexpr auto thermalConductivityMoistureAndTemperatureDependent{false};

    SimulationProperties::Instance().setCalculationParameters(
      excludeWaterLiquidTransportation,
      excludeHeatOfEvaporation,
      excludeCapillaryConduction,
      excludeVaporDiffusionConduction,
      thermalConductivityMoistureAndTemperatureDependent);

    HygroThermFEM::MultiDomain multiDomain;

    const auto & material = multiDomain.materials().createSolidMaterial(TestHelper::CottaerSandstone());

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

    std::vector<std::vector<double>> correctWaterContentSolution{
      {62.292479, 62.292479, 62.999995, 62.999995, 62.292479, 62.292479},
      {61.562200, 61.562200, 62.999994, 62.999994, 61.562200, 61.562200},
      {60.809466, 60.809466, 62.999987, 62.999987, 60.809466, 60.809466},
      {60.033331, 60.033331, 62.999973, 62.999973, 60.033331, 60.033331},
      {59.232351, 59.232351, 62.999953, 62.999953, 59.232351, 59.232351},
      {58.404674, 58.404674, 62.999927, 62.999927, 58.404674, 58.404674},
      {57.547999, 57.547999, 62.999894, 62.999894, 57.547999, 57.547999},
      {56.659498, 56.659498, 62.999853, 62.999853, 56.659498, 56.659498},
      {55.735687, 55.735687, 62.999806, 62.999806, 55.735687, 55.735687},
      {54.772254, 54.772254, 62.999751, 62.999751, 54.772254, 54.772254},
      {53.763793, 53.763793, 62.999689, 62.999689, 53.763793, 53.763793},
      {52.703419, 52.703419, 62.999618, 62.999618, 52.703419, 52.703419},
      {51.582160, 51.582160, 62.999539, 62.999539, 51.582160, 51.582160},
      {50.387978, 50.387978, 62.999450, 62.999450, 50.387978, 50.387978},
      {49.104084, 49.104084, 62.999351, 62.999351, 49.104084, 49.104084},
      {47.705807, 47.705807, 62.999241, 62.999241, 47.705807, 47.705807},
      {46.151972, 46.151972, 62.999118, 62.999118, 46.151972, 46.151972},
      {44.373744, 44.373744, 62.998980, 62.998980, 44.373744, 44.373744},
      {42.509151, 42.509151, 62.998827, 62.998827, 42.509151, 42.509151},
      {40.604655, 40.604655, 62.998660, 62.998660, 40.604655, 40.604655},
      {38.657527, 38.657527, 62.998477, 62.998477, 38.657527, 38.657527},
      {36.664606, 36.664606, 62.998279, 62.998279, 36.664606, 36.664606},
      {34.622323, 34.622323, 62.998065, 62.998065, 34.622323, 34.622323},
      {32.526621, 32.526621, 62.997834, 62.997834, 32.526621, 32.526621}};

    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{
      {0.114108, 0.114108, 0.078490, 0.078490, 0.114108, 0.114108},
      {0.188522, 0.188522, 0.154179, 0.154179, 0.188522, 0.188522},
      {0.245899, 0.245899, 0.217273, 0.217273, 0.245899, 0.245899},
      {0.292055, 0.292055, 0.268715, 0.268715, 0.292055, 0.292055},
      {0.329682, 0.329682, 0.310652, 0.310652, 0.329682, 0.329682},
      {0.360619, 0.360619, 0.345020, 0.345020, 0.360619, 0.360619},
      {0.386289, 0.386289, 0.373401, 0.373401, 0.386289, 0.386289},
      {0.407819, 0.407819, 0.397067, 0.397067, 0.407819, 0.407819},
      {0.426112, 0.426112, 0.417035, 0.417035, 0.426112, 0.426112},
      {0.441892, 0.441892, 0.434122, 0.434122, 0.441892, 0.441892},
      {0.455743, 0.455743, 0.448982, 0.448982, 0.455743, 0.455743},
      {0.468138, 0.468138, 0.462148, 0.462148, 0.468138, 0.468138},
      {0.479468, 0.479468, 0.474054, 0.474054, 0.479468, 0.479468},
      {0.490061, 0.490061, 0.485062, 0.485062, 0.490061, 0.490061},
      {0.500202, 0.500202, 0.495481, 0.495481, 0.500202, 0.500202},
      {0.510155, 0.510155, 0.505590, 0.505590, 0.510155, 0.510155},
      {0.520196, 0.520196, 0.515667, 0.515667, 0.520196, 0.520196},
      {0.530647, 0.530647, 0.526021, 0.526021, 0.530647, 0.530647},
      {0.541490, 0.541490, 0.536727, 0.536727, 0.541490, 0.541490},
      {0.552679, 0.552679, 0.547780, 0.547780, 0.552679, 0.552679},
      {0.564207, 0.564207, 0.559177, 0.559177, 0.564207, 0.564207},
      {0.576075, 0.576075, 0.570917, 0.570917, 0.576075, 0.576075},
      {0.588291, 0.588291, 0.583005, 0.583005, 0.588291, 0.588291},
      {0.600869, 0.600869, 0.595454, 0.595454, 0.600869, 0.600869}};

    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-6);
}
