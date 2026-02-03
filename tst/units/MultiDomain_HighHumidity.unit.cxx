#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::NodePool;
using HygroThermFEM::MaterialPool;

class ObserveSimulationProgrees : public Timesteps::TimestepObserver
{
public:
    void levelChanged(unsigned divisionLevel, unsigned) override
    {
        // No need to notify simulation at level zero
        if(divisionLevel > 0)
        {
            m_SimulationCalls.at(divisionLevel) += 1;
        }
    }

    [[nodiscard]] unsigned getLevelOne() const {return m_SimulationCalls.at(1);}
    [[nodiscard]] unsigned getLevelTwo() const {return m_SimulationCalls.at(2);}
    [[nodiscard]] unsigned getLevelThree() const {return m_SimulationCalls.at(3);}

private:
    // Map will simply keep track of how many times simulation was called
    // at given division level
    std::map<unsigned, unsigned> m_SimulationCalls{{1,0}, {2,0}, {3,0}};
};

class MultiDomain_HighHumidity : public testing::Test
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

TEST_F(MultiDomain_HighHumidity, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Simple two elements example with moisture transfer.");

    std::vector<double> gridXCoordinates{0.15, 0.05, 0.00};

    const auto domainTemperature = 0.0;
    const auto domainHumidity = 0.999;
    const auto domainPressure = 101325.0;
    const auto liquidPercent = 1.0;

    HygroThermFEM::State state(domainTemperature, domainHumidity, domainPressure, liquidPercent);
    size_t nodeIndex = 0;
    for(auto val : gridXCoordinates)
    {
        ++nodeIndex;
        NodePool::Instance().createNode(nodeIndex, val, 0.05, state);
        ++nodeIndex;
        NodePool::Instance().createNode(nodeIndex, val, 0.00, state);
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

    ObserveSimulationProgrees progressThermal;
    domain.subscribeThermal(&progressThermal);

    ObserveSimulationProgrees progressMoisture;
    domain.subscribeMoisture(&progressMoisture);

    /// Create elements
    for(size_t i = 1; i <= (NodePool::Instance().maxIndex() - 2) / 2; ++i)
    {
        const auto node1 = 2u * i + 1u;
        const auto node2 = 2u * i + 2u;
        const auto node3 = 2u * i;
        const auto node4 = 2u * i - 1u;
        domain.createElement(node2, node3, node4, node1, material.name());
    }

    // Create Boundary Conditions
    const auto hc = 10.0;
    const auto airTemperature = 20.0;
    const auto airHumidity = 1.0;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff{airTemperature, hc, airHumidity};

    domain.createBC_FixedHc(5, 6, bcCoeff);

    const auto dTime = 3600;
    const auto nSteps = 2;

    auto temperatures = NodePool::Instance().properties(HygroThermFEM::Variable::temperature);
    auto humidities = NodePool::Instance().properties(HygroThermFEM::Variable::humidity);
    std::vector<std::vector<double>> temperatureSolution;
    std::vector<double> temperatureError;
    std::vector<std::vector<double>> waterContentSolution;
    std::vector<double> humidityError;
    size_t timestepIndex{0};

    for(auto i = 0; i < nSteps; ++i)
    {
        auto aSolution = domain.transient(temperatures, humidities, dTime, timestepIndex);
        temperatureSolution.push_back(aSolution.temperature);
        temperatureError.push_back(aSolution.temperatureError);
        waterContentSolution.push_back(aSolution.waterContent);
        humidityError.push_back(aSolution.humidityError);
        temperatures = aSolution.temperature;
        humidities = aSolution.humidity;
        ++timestepIndex;
    }

    const std::vector<double> correctHumidityError{4.131868e-07, 1.509739e-06};
    const std::vector<std::vector<double>> correctWaterContentSolution{
            {121.994944,121.994944,122.127524,122.127524,122.292494,122.292494},
            {123.768705,123.768705,123.876207,123.876207,124.010072,124.010072}};

    EXPECT_EQ(waterContentSolution.size(), correctWaterContentSolution.size());

    for(auto i = 0u; i < waterContentSolution.size(); ++i)
    {
        EXPECT_NEAR(correctHumidityError[i], humidityError[i], 1e-10);
        for(auto j = 0u; j < waterContentSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctWaterContentSolution[i][j], waterContentSolution[i][j], 1e-6);
        }
    }

    const std::vector<double> correctTemperatureError{0.005612, 0.013515};
    const std::vector<std::vector<double>> correctTemperatureSolution{
            {0.696524,0.696524,2.284785,2.284785,6.987855,6.987855},
            {1.881866,1.881866,4.601701,4.601701,9.950096,9.950096}};

    EXPECT_EQ(temperatureSolution.size(), correctTemperatureSolution.size());

    for(auto i = 0u; i < correctTemperatureSolution.size(); ++i)
    {
        EXPECT_NEAR(correctTemperatureError[i], temperatureError[i], 1e-6);
        for(auto j = 0u; j < correctTemperatureSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctTemperatureSolution[i][j], temperatureSolution[i][j], 1e-6);
        }
    }

    // Checking number of iterations within subiterations

    auto lvlOneMoisture = progressMoisture.getLevelOne();
    EXPECT_EQ(lvlOneMoisture, 7u);

    auto lvlTwoMoisture = progressMoisture.getLevelTwo();
    EXPECT_EQ(lvlTwoMoisture, 7u);

    auto lvlThreeMoisture = progressMoisture.getLevelThree();
    EXPECT_EQ(lvlThreeMoisture, 7007u);

    auto lvlOneThermal = progressThermal.getLevelOne();
    EXPECT_EQ(lvlOneThermal, 0u);

    auto lvlTwoThermal = progressThermal.getLevelTwo();
    EXPECT_EQ(lvlTwoThermal, 0u);

    auto lvlThreeThermal = progressThermal.getLevelThree();
    EXPECT_EQ(lvlThreeThermal, 0u);
}
