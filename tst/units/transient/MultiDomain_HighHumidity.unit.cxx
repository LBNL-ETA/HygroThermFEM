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

    const std::vector correctHumidityError{6.576671e-09, 5.892850e-09};
    const std::vector<std::vector<double>> correctWaterContentSolution{
      {55.076174, 55.077641, 47.684512, 47.687252, 38.379533, 38.381033},
      {53.000743, 52.947827, 45.616675, 45.514793, 32.950396, 32.883968}};

    TestHelper::expectNear(correctHumidityError, humidityError, 1e-6);
    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-3);

    const std::vector correctTemperatureError{4.297336e-07, 1.857763e-07};
    const std::vector<std::vector<double>> correctTemperatureSolution{
      {1.058283, 1.058283, 3.461718, 3.461712, 9.567513, 9.567507},
      {3.008612, 3.008642, 6.604594, 6.604754, 12.124579, 12.124799}};

    TestHelper::expectNear(correctTemperatureError, temperatureError, 1e-7);
    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-3);

    // Checking number of iterations within subiterations
    EXPECT_EQ(progressMoisture.calls(), (std::vector<unsigned>{6u, 0u, 0u}));
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

    const std::vector correctHumidityError{4.828904e-09, 2.026291e-06};
    const std::vector<std::vector<double>> correctWaterContentSolution{
      {60.955985, 60.955984, 51.108232, 51.108232, 35.161279, 35.161279},
      {76.041154, 76.041159, 180.000000, 180.000000, 35.577946, 35.577903}};

    TestHelper::expectNear(correctHumidityError, humidityError, 1e-6);
    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-3);

    const std::vector correctTemperatureError{2.530341e-08, 2.013575e-07};
    const std::vector<std::vector<double>> correctTemperatureSolution{
      {29.460719, 29.460719, 27.976497, 27.976497, 24.218553, 24.218553},
      {34.631913, 34.631913, 42.885187, 42.885188, 53.838928, 53.838936}};

    TestHelper::expectNear(correctTemperatureError, temperatureError, 1e-6);
    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-3);

    // Checking number of iterations within subiterations
    EXPECT_EQ(progressMoisture.calls(), (std::vector<unsigned>{5u, 0u, 0u}));
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

    const std::vector<double> correctHumidityError{7.423032e-09, 2.439338e-07};
    const std::vector<std::vector<double>> correctWaterContentSolution{
      {97.071668, 75.241164, 63.629275, 61.322039, 37.362987, 36.859224},
      {0.000000, 5.300000, 0.000000, 5.300000, 2.619785, 2.458302}};

    TestHelper::expectNear(correctHumidityError, humidityError, 1e-6);
    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-3);

    const std::vector<double> correctTemperatureError{6.345433e-08, 7.502573e-08};
    const std::vector<std::vector<double>> correctTemperatureSolution{
      {75.569058, 75.510763, 64.128790, 64.111068, 37.799341, 37.805086},
      {68.739455, 68.801127, 59.512914, 59.543755, 50.117837, 50.289292}};

    TestHelper::expectNear(correctTemperatureError, temperatureError, 1e-6);
    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-3);

    // Checking number of iterations within subiterations
    EXPECT_EQ(progressMoisture.calls(), (std::vector<unsigned>{3u, 0u, 0u}));
    EXPECT_EQ(progressThermal.calls(), (std::vector<unsigned>{0u, 0u, 0u}));
}
