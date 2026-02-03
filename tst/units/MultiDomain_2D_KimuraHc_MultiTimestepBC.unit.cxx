#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::NodePool;
using HygroThermFEM::MaterialPool;

class MultiDomain_2D_KimuraHc_MultiTimestepBC : public testing::Test
{
protected:
    void SetUp() override
    {}

    void TearDown() override
    {
        NodePool::Instance().clear();
        MaterialPool::Instance().clear();
    }
};

TEST_F(MultiDomain_2D_KimuraHc_MultiTimestepBC, TestExample_1)
{
    // Enter nodes. Arguments are: node number, x-coordinate, y-coordinate

    std::vector<double> gridXCoordinates{0, 0.05, 0.1};

    const double initialTemperature = 0.0;
    const double initialMoistureContent = 0.0;
    const double initialPressure = 101325;
    constexpr auto liquidPercent = 1.0;

    auto state = HygroThermFEM::State(
      initialTemperature, initialMoistureContent, initialPressure, liquidPercent);
    size_t nodeIndex = 0;
    for(auto val : gridXCoordinates)
    {
        ++nodeIndex;
        NodePool::Instance().createNode(nodeIndex, val, 0.00, state);
        ++nodeIndex;
        NodePool::Instance().createNode(nodeIndex, val, 0.05, state);
    }

    // Material Properties (Cottaer Sandstone)
    constexpr double thermalConductivityDry{1.8};
    constexpr double density{2050.0};
    constexpr double porosity{0.22};
    constexpr double specificHeatCapacityDry{850.0};
    constexpr double diffusionResistanceFactor{15.0};
    const std::vector<FenestrationCommon::point> thermalConductivityMoistureDependent = {
      {0.0, 1.8}, {180, 1.8}};
    constexpr double thermalConductivityMeasuredAtTemperature{0};
    const std::vector<FenestrationCommon::point> thermalConductivityTemperatureDependent = {
      {0.0, 1.8}, {1, 1.8}};
    constexpr double thermalConductivityMeasuredAtHumidity{0};
    const std::vector<FenestrationCommon::point> liquidTransportationCurve = {{0, 0},
                                                                              {27, 1E-8},
                                                                              {45, 1.1E-8},
                                                                              {90, 2E-8},
                                                                              {126, 3.5E-8},
                                                                              {144, 5E-8},
                                                                              {162, 1E-7},
                                                                              {171, 2E-7},
                                                                              {180, 7E-7}};

    const std::vector<FenestrationCommon::point> moistureStorageFunction = {{0, 0},
                                                                            {0.5, 5.3},
                                                                            {0.65, 8.4},
                                                                            {0.8, 12},
                                                                            {0.93, 17},
                                                                            {0.95, 25},
                                                                            {0.99, 63},
                                                                            {0.995, 83},
                                                                            {0.999, 120},
                                                                            {1, 180}};

    auto & material =
      MaterialPool::Instance().createSolidMaterial("Cottaer Sandstone",
                                                   thermalConductivityDry,
                                                   density,
                                                   porosity,
                                                   specificHeatCapacityDry,
                                                   diffusionResistanceFactor,
                                                   thermalConductivityMoistureDependent,
                                                   thermalConductivityMeasuredAtTemperature,
                                                   thermalConductivityTemperatureDependent,
                                                   thermalConductivityMeasuredAtHumidity,
                                                   liquidTransportationCurve,
                                                   moistureStorageFunction);

    HygroThermFEM::MultiDomain domain;

    /// Create elements
    for(size_t i = 1; i <= (NodePool::Instance().maxIndex() - 2) / 2; ++i)
    {
        const auto node1 = 2u * i + 1u;
        const auto node2 = 2u * i + 2u;
        const auto node3 = 2u * i;
        const auto node4 = 2u * i - 1u;
        domain.createElement(node1, node2, node3, node4, material.name());
    }

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

    domain.createBC_KimuraHc(1, 2, bcCoeff);

    constexpr auto dTime = 3600;
    constexpr auto nSteps = 10;

    auto temperatures = NodePool::Instance().properties(HygroThermFEM::Variable::temperature);
    auto humidities = NodePool::Instance().properties(HygroThermFEM::Variable::humidity);
    std::vector<std::vector<double>> temperatureSolution;
    std::vector<std::vector<double>> waterContentSolution;
    size_t timestepIndex{0u};

    for(auto i = 0; i < nSteps; ++i)
    {
        auto aSolution = domain.transient(temperatures, humidities, dTime, timestepIndex);
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
