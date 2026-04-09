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
    std::vector<double> temperatureError;
    std::vector<std::vector<double>> waterContentSolution;
    std::vector<double> humidityError;

    for(auto i = 0; i < nSteps; ++i)
    {
        auto aSolution = multiDomain.transient(temperatures, humidities, dTime);
        temperatureSolution.push_back(aSolution.temperature);
        temperatureError.push_back(aSolution.temperatureError);
        waterContentSolution.push_back(aSolution.waterContent);
        humidityError.push_back(aSolution.humidityError);
        temperatures = aSolution.temperature;
        humidities = aSolution.humidity;
    }

    const std::vector<std::vector<double>> correctWaterContentSolution{
      {62.535544, 62.535544, 62.763797, 62.763797, 62.535544, 62.535544},
      {62.138115, 62.138115, 62.445593, 62.445593, 62.138115, 62.138115},
      {61.758606, 61.758606, 62.096715, 62.096715, 61.758606, 61.758606},
      {61.381410, 61.381410, 61.731897, 61.731897, 61.381410, 61.381410},
      {61.000119, 61.000119, 61.357414, 61.357414, 61.000119, 61.000119},
      {60.611503, 60.611503, 60.974477, 60.974477, 60.611503, 60.611503},
      {60.216359, 60.216359, 60.583663, 60.583663, 60.216359, 60.216359},
      {59.813802, 59.813802, 60.185403, 60.185403, 59.813802, 59.813802},
      {59.403529, 59.403529, 59.779439, 59.779439, 59.403529, 59.403529},
      {58.985188, 58.985188, 59.365449, 59.365449, 58.985188, 58.985188},
      {58.557797, 58.557797, 58.942903, 58.942903, 58.557797, 58.557797},
      {58.122426, 58.122426, 58.510519, 58.510519, 58.122426, 58.122426},
      {57.676498, 57.676498, 58.070572, 58.070572, 57.676498, 57.676498},
      {57.221385, 57.221385, 57.619007, 57.619007, 57.221385, 57.221385},
      {56.754651, 56.754651, 57.158460, 57.158460, 56.754651, 56.754651},
      {56.277275, 56.277275, 56.685661, 56.685661, 56.277275, 56.277275},
      {55.786413, 55.786413, 56.200782, 56.200782, 55.786413, 55.786413},
      {55.283368, 55.283368, 55.702311, 55.702311, 55.283368, 55.283368},
      {54.764931, 54.764931, 55.190243, 55.190243, 54.764931, 54.764931},
      {54.231117, 54.231117, 54.661896, 54.661896, 54.231117, 54.231117},
      {53.679199, 53.679199, 54.116841, 54.116841, 53.679199, 53.679199},
      {53.108478, 53.108478, 53.552086, 53.552086, 53.108478, 53.108478},
      {52.515204, 52.515204, 52.966922, 52.966922, 52.515204, 52.515204},
      {51.899108, 51.899108, 52.356291, 52.356291, 51.899108, 51.899108}};

    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    const std::vector<std::vector<double>> correctTemperatureSolution{
      {0.113693, 0.113693, 0.078277, 0.078277, 0.113693, 0.113693},
      {0.185509, 0.185509, 0.154952, 0.154952, 0.185509, 0.185509},
      {0.241347, 0.241347, 0.218105, 0.218105, 0.241347, 0.241347},
      {0.286424, 0.286424, 0.268771, 0.268771, 0.286424, 0.286424},
      {0.322921, 0.322921, 0.309350, 0.309350, 0.322921, 0.322921},
      {0.352482, 0.352482, 0.341927, 0.341927, 0.352482, 0.352482},
      {0.376468, 0.376468, 0.368174, 0.368174, 0.376468, 0.376468},
      {0.396000, 0.396000, 0.389419, 0.389419, 0.396000, 0.396000},
      {0.411993, 0.411993, 0.406718, 0.406718, 0.411993, 0.411993},
      {0.425180, 0.425180, 0.420908, 0.420908, 0.425180, 0.425180},
      {0.436154, 0.436154, 0.432654, 0.432654, 0.436154, 0.436154},
      {0.445386, 0.445386, 0.442482, 0.442482, 0.445386, 0.445386},
      {0.453253, 0.453253, 0.450809, 0.450809, 0.453253, 0.453253},
      {0.460055, 0.460055, 0.457966, 0.457966, 0.460055, 0.460055},
      {0.466031, 0.466031, 0.464215, 0.464215, 0.466031, 0.466031},
      {0.471370, 0.471370, 0.469764, 0.469764, 0.471370, 0.471370},
      {0.476226, 0.476226, 0.474780, 0.474780, 0.476226, 0.476226},
      {0.480717, 0.480717, 0.479392, 0.479392, 0.480717, 0.480717},
      {0.484944, 0.484944, 0.483708, 0.483708, 0.484944, 0.484944},
      {0.488983, 0.488983, 0.487811, 0.487811, 0.488983, 0.488983},
      {0.492899, 0.492899, 0.491772, 0.491772, 0.492899, 0.492899},
      {0.496746, 0.496746, 0.495647, 0.495647, 0.496746, 0.496746},
      {0.500570, 0.500570, 0.499485, 0.499485, 0.500570, 0.500570},
      {0.504407, 0.504407, 0.503326, 0.503326, 0.504407, 0.504407}};

    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-6);
}
