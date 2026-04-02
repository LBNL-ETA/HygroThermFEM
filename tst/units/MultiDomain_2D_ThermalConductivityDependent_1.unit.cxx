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

    std::vector<std::vector<double>> correctWaterContentSolution{
        {13.412206,13.412206,25.252290,25.252290,13.412206,13.412206},
        {12.722048,12.722048,25.246485,25.246485,12.722048,12.722048},
        {12.026558,12.026558,25.239667,25.239667,12.026558,12.026558},
        {11.580200,11.580200,25.231835,25.231835,11.580200,11.580200},
        {11.142305,11.142305,25.222985,25.222985,11.142305,11.142305},
        {10.703634,10.703634,25.213108,25.213108,10.703634,10.703634},
        {10.264813,10.264813,25.202190,25.202190,10.264813,10.264813},
        {9.830200,9.830200,25.190232,25.190232,9.830200,9.830200},
        {9.400823,9.400823,25.177234,25.177234,9.400823,9.400823},
        {8.976953,8.976953,25.163191,25.163191,8.976953,8.976953},
        {8.558806,8.558806,25.148102,25.148102,8.558806,8.558806},
        {8.181753,8.181753,25.131962,25.131962,8.181753,8.181753},
        {7.831947,7.831947,25.114770,25.114770,7.831947,7.831947},
        {7.487427,7.487427,25.096522,25.096522,7.487427,7.487427},
        {7.148270,7.148270,25.077215,25.077215,7.148270,7.148270},
        {6.814501,6.814501,25.056845,25.056845,6.814501,6.814501},
        {6.485884,6.485884,25.035408,25.035408,6.485884,6.485884},
        {6.162483,6.162483,25.012900,25.012900,6.162483,6.162483},
        {5.844354,5.844354,24.995502,24.995502,5.844354,5.844354},
        {5.531544,5.531544,24.985118,24.985118,5.531544,5.531544},
        {5.261069,5.261069,24.974279,24.974279,5.261069,5.261069},
        {5.106150,5.106150,24.962983,24.962983,5.106150,5.106150},
        {4.954018,4.954018,24.951228,24.951228,4.954018,4.954018},
        {4.804690,4.804690,24.939015,24.939015,4.804690,4.804690}};

    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{
        {0.387303,0.387303,0.281075,0.281075,0.387303,0.387303},
        {0.644536,0.644536,0.573262,0.573262,0.644536,0.644536},
        {0.875623,0.875623,0.824429,0.824429,0.875623,0.875623},
        {1.083743,1.083743,1.044239,1.044239,1.083743,1.083743},
        {1.272934,1.272934,1.240924,1.240924,1.272934,1.272934},
        {1.447448,1.447448,1.420495,1.420495,1.447448,1.447448},
        {1.610878,1.610878,1.587478,1.587478,1.610878,1.610878},
        {1.765821,1.765821,1.745038,1.745038,1.765821,1.765821},
        {1.914285,1.914285,1.895487,1.895487,1.914285,1.914285},
        {2.057811,2.057811,2.040559,2.040559,2.057811,2.057811},
        {2.197570,2.197570,2.181552,2.181552,2.197570,2.197570},
        {2.334436,2.334436,2.319427,2.319427,2.334436,2.334436},
        {2.469071,2.469071,2.454904,2.454904,2.469071,2.469071},
        {2.601976,2.601976,2.588525,2.588525,2.601976,2.601976},
        {2.733518,2.733518,2.720688,2.720688,2.733518,2.733518},
        {2.863972,2.863972,2.851689,2.851689,2.863972,2.863972},
        {2.993563,2.993563,2.981764,2.981764,2.993563,2.993563},
        {3.122448,3.122448,3.111085,3.111085,3.122448,3.122448},
        {3.250734,3.250734,3.239769,3.239769,3.250734,3.250734},
        {3.378493,3.378493,3.367893,3.367893,3.378493,3.378493},
        {3.505757,3.505757,3.495496,3.495496,3.505757,3.505757},
        {3.632520,3.632520,3.622578,3.622578,3.632520,3.632520},
        {3.758798,3.758798,3.749153,3.749153,3.758798,3.758798},
        {3.884588,3.884588,3.875224,3.875224,3.884588,3.884588}};

    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-6);
}
