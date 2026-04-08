#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "TestMaterials.hxx"
#include "TestHelpers.hxx"

using HygroThermFEM::Nodes;

namespace
{
    //! Counts how many times each non-zero subdivision level fires.
    //! Index 0 corresponds to division level 1, index 1 to level 2, etc.
    class ObserveSimulationProgress : public Timesteps::TimestepObserver
    {
    public:
        void levelChanged(unsigned divisionLevel, unsigned) override
        {
            // Level zero is the initial timestep attempt and is not counted.
            if(divisionLevel == 0)
            {
                return;
            }
            const size_t idx = divisionLevel - 1;
            if(idx >= m_SimulationCalls.size())
            {
                m_SimulationCalls.resize(idx + 1, 0u);
            }
            ++m_SimulationCalls[idx];
        }

        [[nodiscard]] const std::vector<unsigned> & calls() const
        {
            return m_SimulationCalls;
        }

    private:
        std::vector<unsigned> m_SimulationCalls{0u, 0u, 0u};
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
    EXPECT_EQ(progressMoisture.calls(), (std::vector<unsigned>{8u, 7u, 3u}));
    EXPECT_EQ(progressThermal.calls(), (std::vector<unsigned>{0u, 0u, 0u}));
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
    EXPECT_EQ(progressMoisture.calls(), (std::vector<unsigned>{4u, 4u, 1u}));
    EXPECT_EQ(progressThermal.calls(), (std::vector<unsigned>{0u, 0u, 0u}));
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
    EXPECT_EQ(progressMoisture.calls(), (std::vector<unsigned>{25u, 25u, 22u}));
    EXPECT_EQ(progressThermal.calls(), (std::vector<unsigned>{0u, 0u, 0u}));
}
