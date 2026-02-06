#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "SlabCreator.hxx"
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
      {4.630345, 4.630345, 0.007231, 0.007231, 2.1e-05, 2.1e-05},
      {5.851888, 5.851888, 0.019015, 0.019015, 9.8e-05, 9.8e-05},
      {4.915636, 4.915636, 0.032528, 0.032528, 0.000271, 0.000271},
      {3.251349, 3.251349, 0.044297, 0.044297, 0.000574, 0.000574},
      {1.761459, 1.761459, 0.052298, 0.052298, 0.001029, 0.001029},
      {1.320202, 1.320202, 0.058951, 0.058951, 0.001628, 0.001628},
      {1.0834, 1.0834, 0.064677, 0.064677, 0.002349, 0.002349},
      {0.940131, 0.940131, 0.069671, 0.069671, 0.00316, 0.00316},
      {0.847121, 0.847121, 0.074042, 0.074042, 0.004026, 0.004026},
      {0.784377, 0.784377, 0.077867, 0.077867, 0.004913, 0.004913}};

    EXPECT_EQ(waterContentSolution.size(), correctWaterContentSolution.size());

    for(auto i = 0u; i < correctWaterContentSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctWaterContentSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctWaterContentSolution[i][j], waterContentSolution[i][j], 1e-6);
        }
    }

    std::vector<std::vector<double>> correctTemperatureSolution{
      {11.525284, 11.525284, 5.985203, 5.985203, 4.479496, 4.479496},
      {16.118553, 16.118553, 11.048928, 11.048928, 9.396227, 9.396227},
      {19.439666, 19.439666, 15.188898, 15.188898, 13.731587, 13.731587},
      {23.656529, 23.656529, 19.395658, 19.395658, 17.970676, 17.970676},
      {27.908635, 27.908635, 23.631663, 23.631663, 22.207435, 22.207435},
      {28.64284, 28.64284, 26.048609, 26.048609, 25.082191, 25.082191},
      {29.089923, 29.089923, 27.501981, 27.501981, 26.893145, 26.893145},
      {28.97192, 28.97192, 28.185653, 28.185653, 27.860414, 27.860414},
      {28.33519, 28.33519, 28.220333, 28.220333, 28.129715, 28.129715},
      {27.262137, 27.262137, 27.710112, 27.710112, 27.81561, 27.81561}};

    EXPECT_EQ(temperatureSolution.size(), correctTemperatureSolution.size());

    for(auto i = 0u; i < correctTemperatureSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctTemperatureSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctTemperatureSolution[i][j], temperatureSolution[i][j], 1e-6);
        }
    }
}
