#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "SlabCreator.hxx"
#include "TestMaterials.hxx"

TEST(MultiDomain_2D_YazdanianKlemsHc_MultiTimestepBC, TestExample_1)
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

    // Variable boundary conditions (temperature, humidity, wind speed and wind direction) over ten timesteps.
    const std::vector<HygroThermFEM::YazdanianKlemsCoefficients> bcCoeff{
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

    multiDomain.createBC_YazdanianKlemsHc(1, 2, bcCoeff);

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
      {4.119317, 4.119317, 0.005434, 0.005434, 1.4e-05, 1.4e-05},
      {6.250393, 6.250393, 0.014202, 0.014202, 5.7e-05, 5.7e-05},
      {7.302342, 7.302342, 0.02489, 0.02489, 0.000144, 0.000144},
      {6.784836, 6.784836, 0.036153, 0.036153, 0.000285, 0.000285},
      {5.107654, 5.107654, 0.046459, 0.046459, 0.000479, 0.000479},
      {4.352747, 4.352747, 0.055702, 0.055702, 0.000724, 0.000724},
      {3.673172, 3.673172, 0.063756, 0.063756, 0.001013, 0.001013},
      {3.107442, 3.107442, 0.070658, 0.070658, 0.001337, 0.001337},
      {2.663355, 2.663355, 0.076545, 0.076545, 0.001687, 0.001687},
      {2.333544, 2.333544, 0.081553, 0.081553, 0.002055, 0.002055}};

    EXPECT_EQ(waterContentSolution.size(), correctWaterContentSolution.size());

    for(auto i = 0u; i < correctWaterContentSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctWaterContentSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctWaterContentSolution[i][j], waterContentSolution[i][j], 1e-6);
        }
    }

    std::vector<std::vector<double>> correctTemperatureSolution{
      {6.290785, 6.290785, 3.267036, 3.267036, 2.445143, 2.445143},
      {8.607262, 8.607262, 5.931511, 5.931511, 5.054433, 5.054433},
      {9.762976, 9.762976, 7.804778, 7.804778, 7.112854, 7.112854},
      {10.832759, 10.832759, 9.285138, 9.285138, 8.738631, 8.738631},
      {11.490438, 11.490438, 10.35778, 10.35778, 9.950421, 9.950421},
      {11.749038, 11.749038, 11.026077, 11.026077, 10.755444, 10.755444},
      {11.894196, 11.894196, 11.440727, 11.440727, 11.2683, 11.2683},
      {11.837921, 11.837921, 11.623755, 11.623755, 11.534303, 11.534303},
      {11.567187, 11.567187, 11.582029, 11.582029, 11.569995, 11.569995},
      {11.075398, 11.075398, 11.316698, 11.316698, 11.380395, 11.380395}};

    EXPECT_EQ(temperatureSolution.size(), correctTemperatureSolution.size());

    for(auto i = 0u; i < correctTemperatureSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctTemperatureSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctTemperatureSolution[i][j], temperatureSolution[i][j], 1e-6);
        }
    }
}
