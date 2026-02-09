#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "SlabCreator.hxx"
#include "TestMaterials.hxx"

using HygroThermFEM::Nodes;
using HygroThermFEM::SimulationProperties;

class MultiDomain_2D_NoMoistureStorageFunction : public testing::Test
{
protected:
    void SetUp() override
    {}

    void TearDown() override
    {
        SimulationProperties::Instance().reset();
    }
};

TEST_F(MultiDomain_2D_NoMoistureStorageFunction, TestExample_1)
{
    constexpr auto excludeWaterLiquidTransportation{false};
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

    // Material Properties to represent example without moisture storage function and liquid
    // transportation
    auto params = TestHelper::CottaerSandstone();
    params.name = "No liquid or vapor transport curves";
    // Zero out moisture storage and liquid transport. The sorption curve returns
    // zero water content at all humidities, so max water content = 0.
    // The material validation requires that the thermal conductivity moisture-
    // dependent curve and liquid transport curve have their max water content
    // (x-axis) equal to the sorption curve's max water content (y-axis).
    // Single-point curves get auto-expanded by TabularFunction1D (adds x+1),
    // so we use explicit two-point curves to prevent expansion and keep all
    // max water content values consistently at zero.
    params.sorptionCurve = {{0, 0}, {1, 0}};
    params.liquidTransportCurve = {{0, 0}, {0, 0}};
    params.thermalConductivityMoistureDependent = {{0, 1.8}, {0, 1.8}};
    const auto & material = multiDomain.materials().createSolidMaterial(params);

    constexpr HygroThermFEM::State state({
        .temperature = 0.0,
        .humidity = 0.2,
        .pressure = 101325.0
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
    constexpr auto nSteps = 4;

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

    std::vector<std::vector<double>> correctWaterContentSolution{{0, 0, 0, 0, 0, 0},
                                                                 {0, 0, 0, 0, 0, 0},
                                                                 {0, 0, 0, 0, 0, 0},
                                                                 {0, 0, 0, 0, 0, 0}};

    EXPECT_EQ(waterContentSolution.size(), correctWaterContentSolution.size());

    for(auto i = 0u; i < waterContentSolution.size(); ++i)
    {
        for(auto j = 0u; j < waterContentSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctWaterContentSolution[i][j], waterContentSolution[i][j], 1e-6);
        }
    }

    std::vector<std::vector<double>> correctTemperatureSolution{
      {1.911518, 1.911518, 1.430638, 1.430638, 1.911518, 1.911518},
      {3.291149, 3.291149, 2.823100, 2.823100, 3.291149, 3.291149},
      {4.411504, 4.411504, 4.011909, 4.011909, 4.411504, 4.411504},
      {5.341402, 5.341402, 5.006941, 5.006941, 5.341402, 5.341402}};

    EXPECT_EQ(temperatureSolution.size(), correctTemperatureSolution.size());

    for(auto i = 0u; i < temperatureSolution.size(); ++i)
    {
        for(auto j = 0u; j < temperatureSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctTemperatureSolution[i][j], temperatureSolution[i][j], 1e-6);
        }
    }
}
