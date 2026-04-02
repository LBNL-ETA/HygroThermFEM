#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "TestMaterials.hxx"
#include "TestHelpers.hxx"

using HygroThermFEM::Nodes;

namespace
{
    class ObserveSimulationProgress : public Timesteps::TimestepObserver
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

        [[nodiscard]] unsigned getLevelOne() const
        {
            return m_SimulationCalls.at(1);
        }
        [[nodiscard]] unsigned getLevelTwo() const
        {
            return m_SimulationCalls.at(2);
        }
        [[nodiscard]] unsigned getLevelThree() const
        {
            return m_SimulationCalls.at(3);
        }

    private:
        // Map will simply keep track of how many times simulation was called
        // at given division level
        std::map<unsigned, unsigned> m_SimulationCalls{{1, 0}, {2, 0}, {3, 0}};
    };
}   // namespace

TEST(MultiDomain_HighHumidity, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Simple two elements example with moisture transfer.");

    HygroThermFEM::MultiDomain multiDomain;

    std::vector<double> gridXCoordinates{0.15, 0.05, 0.00};

    constexpr HygroThermFEM::State state(
      {.temperature = 0.0, .humidity = 0.999, .pressure = 101325.0, .liquidPercent = 1.0});

    for(auto val : gridXCoordinates)
    {
        multiDomain.nodes().createNode({.x = val, .y = 0.05, .state = state});
        multiDomain.nodes().createNode({.x = val, .y = 0.00, .state = state});
    }

    const auto & material =
      multiDomain.materials().createSolidMaterial(TestHelper::CottaerSandstone());

    ObserveSimulationProgress progressThermal;
    multiDomain.subscribeThermal(&progressThermal);

    ObserveSimulationProgress progressMoisture;
    multiDomain.subscribeMoisture(&progressMoisture);

    /// Create elements
    for(size_t i = 1; i <= (multiDomain.nodes().maxIndex() - 2) / 2; ++i)
    {
        const auto node1 = 2u * i + 1u;
        const auto node2 = 2u * i + 2u;
        const auto node3 = 2u * i;
        const auto node4 = 2u * i - 1u;
        multiDomain.createElement({.node1 = node2,
                                   .node2 = node3,
                                   .node3 = node4,
                                   .node4 = node1,
                                   .material = material.name()});
    }

    // Create Boundary Conditions
    constexpr auto hc = 10.0;
    constexpr auto airTemperature = 20.0;
    constexpr auto airHumidity = 1.0;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff{airTemperature, hc, airHumidity};

    multiDomain.createBC_FixedHc(5, 6, bcCoeff);

    constexpr auto dTime = 3600;
    constexpr auto nSteps = 2;

    auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
    auto humidities = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);
    std::vector<std::vector<double>> temperatureSolution;
    std::vector<double> temperatureError;
    std::vector<std::vector<double>> waterContentSolution;
    std::vector<double> humidityError;
    size_t timestepIndex{0};

    for(auto i = 0; i < nSteps; ++i)
    {
        auto aSolution = multiDomain.transient(temperatures, humidities, dTime, timestepIndex);
        temperatureSolution.push_back(aSolution.temperature);
        temperatureError.push_back(aSolution.temperatureError);
        waterContentSolution.push_back(aSolution.waterContent);
        humidityError.push_back(aSolution.humidityError);
        temperatures = aSolution.temperature;
        humidities = aSolution.humidity;
        ++timestepIndex;
    }

    const std::vector correctHumidityError{1.754033e-06, 0.0};
    const std::vector<std::vector<double>> correctWaterContentSolution{
      {180.000000, 180.000000, 180.000000, 180.000000, 27.758519, 27.790584},
      {180.000000, 180.000000, 180.000000, 180.000000, 180.000000, 180.000000}};

    TestHelper::expectNear(correctHumidityError, humidityError, 1e-6);
    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-3);

    const std::vector correctTemperatureError{3.990158e-07, 2.124443e-07};
    const std::vector<std::vector<double>> correctTemperatureSolution{
      {-0.108497, -0.108597, -0.473111, -0.472734, 6.802888, 6.802827},
      {0.535021, 0.534521, 2.337413, 2.339001, 9.005659, 9.005527}};

    TestHelper::expectNear(correctTemperatureError, temperatureError, 1e-7);
    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-3);

    // Checking number of iterations within subiterations

    auto lvlOneMoisture = progressMoisture.getLevelOne();
    EXPECT_EQ(lvlOneMoisture, 8u);

    auto lvlTwoMoisture = progressMoisture.getLevelTwo();
    EXPECT_EQ(lvlTwoMoisture, 7u);

    auto lvlThreeMoisture = progressMoisture.getLevelThree();
    EXPECT_EQ(lvlThreeMoisture, 3u);

    auto lvlOneThermal = progressThermal.getLevelOne();
    EXPECT_EQ(lvlOneThermal, 0u);

    auto lvlTwoThermal = progressThermal.getLevelTwo();
    EXPECT_EQ(lvlTwoThermal, 0u);

    auto lvlThreeThermal = progressThermal.getLevelThree();
    EXPECT_EQ(lvlThreeThermal, 0u);
}

TEST(MultiDomain_HighHumidity, HighHumidityAndTemperature)
{
    SCOPED_TRACE("Begin Test: Simple two elements example with moisture transfer.");

    HygroThermFEM::MultiDomain multiDomain;

    const std::vector gridXCoordinates{0.15, 0.05, 0.00};

    constexpr HygroThermFEM::State state({
        .temperature = 30.0,
        .humidity = 0.9999,
        .pressure = 101325.0,
        .liquidPercent = 1.0
    });

    for(auto val : gridXCoordinates)
    {
        multiDomain.nodes().createNode({.x = val, .y = 0.05, .state = state});
        multiDomain.nodes().createNode({.x = val, .y = 0.00, .state = state});
    }

    const auto & material =
      multiDomain.materials().createSolidMaterial(TestHelper::CottaerSandstone());

    ObserveSimulationProgress progressThermal;
    multiDomain.subscribeThermal(&progressThermal);

    ObserveSimulationProgress progressMoisture;
    multiDomain.subscribeMoisture(&progressMoisture);

    /// Create elements
    for(size_t i = 1; i <= (multiDomain.nodes().maxIndex() - 2) / 2; ++i)
    {
        const auto node1 = 2u * i + 1u;
        const auto node2 = 2u * i + 2u;
        const auto node3 = 2u * i;
        const auto node4 = 2u * i - 1u;
        multiDomain.createElement({.node1 = node2,
                                   .node2 = node3,
                                   .node3 = node4,
                                   .node4 = node1,
                                   .material = material.name()});
    }

    // Create Boundary Conditions
    constexpr auto hc = 10.0;
    constexpr auto airTemperature = 20.0;
    constexpr auto airHumidity = 1.0;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff{airTemperature, hc, airHumidity};

    multiDomain.createBC_FixedHc(5, 6, bcCoeff);

    constexpr auto dTime = 3600;
    constexpr auto nSteps = 2;

    auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
    auto humidities = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);
    std::vector<std::vector<double>> temperatureSolution;
    std::vector<double> temperatureError;
    std::vector<std::vector<double>> waterContentSolution;
    std::vector<double> humidityError;
    size_t timestepIndex{0};

    for(auto i = 0; i < nSteps; ++i)
    {
        auto aSolution = multiDomain.transient(temperatures, humidities, dTime, timestepIndex);
        temperatureSolution.push_back(aSolution.temperature);
        temperatureError.push_back(aSolution.temperatureError);
        waterContentSolution.push_back(aSolution.waterContent);
        humidityError.push_back(aSolution.humidityError);
        temperatures = aSolution.temperature;
        humidities = aSolution.humidity;
        ++timestepIndex;
    }

    const std::vector correctHumidityError{2.657823e-07, 1.591037e-07};
    const std::vector<std::vector<double>> correctWaterContentSolution{
      {81.001587, 81.001587, 180.000000, 180.000000, 180.000000, 180.000000},
      {113.376865, 113.376865, 79.008819, 79.008819, 180.000000, 180.000000}};

    TestHelper::expectNear(correctHumidityError, humidityError, 1e-7);
    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-3);

    const std::vector correctTemperatureError{6.284313e-09, 5.984128e-09};
    const std::vector<std::vector<double>> correctTemperatureSolution{
      {30.985534, 30.985534, 29.173305, 29.173305, 24.758822, 24.758822},
      {30.716207, 30.716207, 30.078265, 30.078265, 27.168913, 27.168913}};

    TestHelper::expectNear(correctTemperatureError, temperatureError, 1e-9);
    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-3);

    // Checking number of iterations within subiterations

    auto lvlOneMoisture = progressMoisture.getLevelOne();
    EXPECT_EQ(lvlOneMoisture, 4u);

    auto lvlTwoMoisture = progressMoisture.getLevelTwo();
    EXPECT_EQ(lvlTwoMoisture, 4u);

    auto lvlThreeMoisture = progressMoisture.getLevelThree();
    EXPECT_EQ(lvlThreeMoisture, 1u);

    auto lvlOneThermal = progressThermal.getLevelOne();
    EXPECT_EQ(lvlOneThermal, 0u);

    auto lvlTwoThermal = progressThermal.getLevelTwo();
    EXPECT_EQ(lvlTwoThermal, 0u);

    auto lvlThreeThermal = progressThermal.getLevelThree();
    EXPECT_EQ(lvlThreeThermal, 0u);
}

TEST(MultiDomain_HighHumidity, ExtremeHumidityAndTemperature)
{
    SCOPED_TRACE("Begin Test: Simple two elements example with moisture transfer.");

    HygroThermFEM::MultiDomain multiDomain;

    std::vector<double> gridXCoordinates{0.15, 0.05, 0.00};

    constexpr HygroThermFEM::State state({
        .temperature = 80.0,
        .humidity = 0.9999,
        .pressure = 101325.0,
        .liquidPercent = 1.0
    });

    for(auto val : gridXCoordinates)
    {
        multiDomain.nodes().createNode({.x = val, .y = 0.05, .state = state});
        multiDomain.nodes().createNode({.x = val, .y = 0.00, .state = state});
    }

    const auto & material =
      multiDomain.materials().createSolidMaterial(TestHelper::CottaerSandstone());

    ObserveSimulationProgress progressThermal;
    multiDomain.subscribeThermal(&progressThermal);

    ObserveSimulationProgress progressMoisture;
    multiDomain.subscribeMoisture(&progressMoisture);

    /// Create elements
    for(size_t i = 1; i <= (multiDomain.nodes().maxIndex() - 2) / 2; ++i)
    {
        const auto node1 = 2u * i + 1u;
        const auto node2 = 2u * i + 2u;
        const auto node3 = 2u * i;
        const auto node4 = 2u * i - 1u;
        multiDomain.createElement({.node1 = node2,
                                   .node2 = node3,
                                   .node3 = node4,
                                   .node4 = node1,
                                   .material = material.name()});
    }

    // Create Boundary Conditions
    constexpr auto hc = 10.0;
    constexpr auto airTemperature = 20.0;
    constexpr auto airHumidity = 1.0;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff{airTemperature, hc, airHumidity};

    multiDomain.createBC_FixedHc(5, 6, bcCoeff);

    constexpr auto dTime = 3600;
    constexpr auto nSteps = 2;

    auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
    auto humidities = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);
    std::vector<std::vector<double>> temperatureSolution;
    std::vector<double> temperatureError;
    std::vector<std::vector<double>> waterContentSolution;
    std::vector<double> humidityError;
    size_t timestepIndex{0};

    for(auto i = 0; i < nSteps; ++i)
    {
        auto aSolution = multiDomain.transient(temperatures, humidities, dTime, timestepIndex);
        temperatureSolution.push_back(aSolution.temperature);
        temperatureError.push_back(aSolution.temperatureError);
        waterContentSolution.push_back(aSolution.waterContent);
        humidityError.push_back(aSolution.humidityError);
        temperatures = aSolution.temperature;
        humidities = aSolution.humidity;
        ++timestepIndex;
    }

    const std::vector<double> correctHumidityError{1.217342e-08, 1.136822e-10};
    const std::vector<std::vector<double>> correctWaterContentSolution{
      {84.901602, 84.901602, 80.217109, 80.217109, 53.239121, 53.239121},
      {180.000000, 180.000000, 180.000000, 180.000000, 180.000000, 180.000000}};

    TestHelper::expectNear(correctHumidityError, humidityError, 1e-8);
    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-3);

    const std::vector<double> correctTemperatureError{1.582014e-07, 5.347030e-08};
    const std::vector<std::vector<double>> correctTemperatureSolution{
      {71.642135, 71.642135, 54.219885, 54.219885, 6.245502, 6.245502},
      {65.864202, 65.864202, 49.926034, 49.926034, 33.073654, 33.073654}};

    TestHelper::expectNear(correctTemperatureError, temperatureError, 1e-7);
    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-3);

    // Checking number of iterations within subiterations

    auto lvlOneMoisture = progressMoisture.getLevelOne();
    EXPECT_EQ(lvlOneMoisture, 25u);

    auto lvlTwoMoisture = progressMoisture.getLevelTwo();
    EXPECT_EQ(lvlTwoMoisture, 25u);

    auto lvlThreeMoisture = progressMoisture.getLevelThree();
    EXPECT_EQ(lvlThreeMoisture, 22u);

    auto lvlOneThermal = progressThermal.getLevelOne();
    EXPECT_EQ(lvlOneThermal, 0u);

    auto lvlTwoThermal = progressThermal.getLevelTwo();
    EXPECT_EQ(lvlTwoThermal, 0u);

    auto lvlThreeThermal = progressThermal.getLevelThree();
    EXPECT_EQ(lvlThreeThermal, 0u);
}
