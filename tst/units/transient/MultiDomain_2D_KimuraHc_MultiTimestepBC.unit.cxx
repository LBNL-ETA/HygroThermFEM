#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "SlabCreator.hxx"
#include "TestHelpers.hxx"
#include "TestMaterials.hxx"

TEST(MultiDomain_2D_KimuraHc_MultiTimestepBC, TestExample_1)
{
    HygroThermFEM::MultiDomain multiDomain;

    const auto & material = multiDomain.materials().createSolidMaterial(TestHelper::CottaerSandstone());

    constexpr HygroThermFEM::State state({
        .temperature = 0.0,
        .humidity = 0.0,
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

    using HygroThermFEM::WindDirection;

    // Variable boundary conditions (temperature, wind speed and humidity) over ten timesteps.
    const std::vector<HygroThermFEM::KimuraCoefficients> bcCoeff{
      {20.0, 0.6, 3.0, WindDirection::Windward},
      {20.0, 0.5, 3.0, WindDirection::Windward},
      {20.0, 0.4, 3.0, WindDirection::Windward},
      {20.0, 0.3, 4.0, WindDirection::Windward},
      {20.0, 0.2, 4.2, WindDirection::Windward},
      {18.0, 0.2, 4.6, WindDirection::Leeward},
      {16.0, 0.2, 5.0, WindDirection::Leeward},
      {14.0, 0.2, 5.3, WindDirection::Leeward},
      {12.0, 0.2, 5.5, WindDirection::Leeward},
      {10.0, 0.2, 5.9, WindDirection::Leeward}};

    multiDomain.createBC_KimuraHc(1, 2, bcCoeff);

    constexpr auto dTime = 3600;
    constexpr auto nSteps = 10;

    auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
    auto humidities = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);
    std::vector<std::vector<double>> temperatureSolution;
    std::vector<std::vector<double>> waterContentSolution;
    size_t timestepIndex{0u};

    for(auto i = 0; i < nSteps; ++i)
    {
        auto aSolution = multiDomain.transient(temperatures, humidities, dTime, timestepIndex);
        temperatureSolution.push_back(aSolution.temperature);
        waterContentSolution.push_back(aSolution.waterContent);
        temperatures = aSolution.temperature;
        humidities = aSolution.humidity;
        ++timestepIndex;
    }

    std::vector<std::vector<double>> correctWaterContentSolution{
      {5.669634, 5.669634, 0.005630, 0.005630, 0.000012, 0.000012},
      {7.785570, 7.785570, 0.015850, 0.015850, 0.000060, 0.000060},
      {6.408183, 6.408183, 0.028616, 0.028616, 0.000181, 0.000181},
      {4.212832, 4.212832, 0.040669, 0.040669, 0.000404, 0.000404},
      {2.342084, 2.342084, 0.049310, 0.049310, 0.000748, 0.000748},
      {1.521257, 1.521257, 0.056198, 0.056198, 0.001243, 0.001243},
      {1.160185, 1.160185, 0.061963, 0.061963, 0.001875, 0.001875},
      {0.957265, 0.957265, 0.066918, 0.066918, 0.002621, 0.002621},
      {0.831072, 0.831072, 0.071219, 0.071219, 0.003448, 0.003448},
      {0.747744, 0.747744, 0.074961, 0.074961, 0.004322, 0.004322}};

    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{
      {12.576657, 12.576657, 6.530750, 6.530750, 4.887802, 4.887802},
      {16.850083, 16.850083, 11.671686, 11.671686, 9.965041, 9.965041},
      {19.839640, 19.839640, 15.687942, 15.687942, 14.248192, 14.248192},
      {23.848185, 23.848185, 19.736655, 19.736655, 18.355863, 18.355863},
      {27.972206, 27.972206, 23.833854, 23.833854, 22.455677, 22.455677},
      {28.594821, 28.594821, 26.126634, 26.126634, 25.203051, 25.203051},
      {29.062018, 29.062018, 27.530512, 27.530512, 26.944912, 26.944912},
      {28.952459, 28.952459, 28.192265, 28.192265, 27.878392, 27.878392},
      {28.317421, 28.317421, 28.215750, 28.215750, 28.130810, 28.130810},
      {27.243859, 27.243859, 27.699143, 27.699143, 27.807677, 27.807677}};

    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-6);
}
