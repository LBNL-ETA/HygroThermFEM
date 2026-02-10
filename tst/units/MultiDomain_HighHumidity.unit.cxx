#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "TestMaterials.hxx"

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

    const std::vector<double> correctHumidityError{2.201043e-09, 6.142575e-09};
    const std::vector<std::vector<double>> correctWaterContentSolution{
      {103.176677, 103.176677, 103.258753, 103.258753, 103.361557, 103.361557},
      {92.753855, 92.753855, 92.595613, 92.595613, 93.380920, 93.380920}};

    EXPECT_EQ(waterContentSolution.size(), correctWaterContentSolution.size());

    for(auto i = 0u; i < waterContentSolution.size(); ++i)
    {
        EXPECT_NEAR(correctHumidityError[i], humidityError[i], 1e-10);
        for(auto j = 0u; j < waterContentSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctWaterContentSolution[i][j], waterContentSolution[i][j], 1e-6);
        }
    }

    const std::vector<double> correctTemperatureError{7.739460e-06, 8.382848e-06};
    const std::vector<std::vector<double>> correctTemperatureSolution{
      {0.783481, 0.783481, 2.452037, 2.452037, 7.203850, 7.203850},
      {2.158483, 2.158483, 4.973106, 4.973106, 10.249844, 10.249844}};

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
    EXPECT_EQ(lvlOneMoisture, 19u);

    auto lvlTwoMoisture = progressMoisture.getLevelTwo();
    EXPECT_EQ(lvlTwoMoisture, 1570u);

    auto lvlThreeMoisture = progressMoisture.getLevelThree();
    EXPECT_EQ(lvlThreeMoisture, 0u);

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

    std::vector<double> gridXCoordinates{0.15, 0.05, 0.00};

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

    const std::vector<double> correctHumidityError{1.210143e-08, 4.116013e-08};
    const std::vector<std::vector<double>> correctWaterContentSolution{
      {105.894452, 105.894452, 105.840409, 105.840409, 105.772928, 105.772928},
      {92.789315, 92.789315, 107.686387, 107.686387, 180.0, 180.0}};

    EXPECT_EQ(waterContentSolution.size(), correctWaterContentSolution.size());

    for(auto i = 0u; i < waterContentSolution.size(); ++i)
    {
        EXPECT_NEAR(correctHumidityError[i], humidityError[i], 1e-10);
        for(auto j = 0u; j < waterContentSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctWaterContentSolution[i][j], waterContentSolution[i][j], 1e-6);
        }
    }

    const std::vector<double> correctTemperatureError{7.080930e-06, 6.976242e-06};
    const std::vector<std::vector<double>> correctTemperatureSolution{
      {29.543444, 29.543444, 28.561120, 28.561120, 25.748466, 25.748466},
      {29.041973, 29.041973, 27.993035, 27.993035, 25.123323, 25.123323}};

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
    EXPECT_EQ(lvlOneMoisture, 19u);

    auto lvlTwoMoisture = progressMoisture.getLevelTwo();
    EXPECT_EQ(lvlTwoMoisture, 1350u);

    auto lvlThreeMoisture = progressMoisture.getLevelThree();
    EXPECT_EQ(lvlThreeMoisture, 0u);

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

    const std::vector<double> correctHumidityError{1.315269e-02, 9.742288e-03};
    const std::vector<std::vector<double>> correctWaterContentSolution{
      {180.0, 180.0, 180.0, 180.0, 180.0, 180.0},
      {180.0, 180.0, 180.0, 180.0, 180.0, 180.0}};

    EXPECT_EQ(waterContentSolution.size(), correctWaterContentSolution.size());

    for(auto i = 0u; i < waterContentSolution.size(); ++i)
    {
        EXPECT_NEAR(correctHumidityError[i], humidityError[i], 1e-4);
        for(auto j = 0u; j < waterContentSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctWaterContentSolution[i][j], waterContentSolution[i][j], 1e-6);
        }
    }

    const std::vector<double> correctTemperatureError{3.802498e-01, 3.575668e-01};
    const std::vector<std::vector<double>> correctTemperatureSolution{
      {79.623285, 79.623285, 78.584134, 78.584134, 75.135358, 75.135358},
      {78.875799, 78.875799, 76.813892, 76.813892, 72.120589, 72.120589}};

    EXPECT_EQ(temperatureSolution.size(), correctTemperatureSolution.size());

    for(auto i = 0u; i < correctTemperatureSolution.size(); ++i)
    {
        EXPECT_NEAR(correctTemperatureError[i], temperatureError[i], 1e-1);
        for(auto j = 0u; j < correctTemperatureSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctTemperatureSolution[i][j], temperatureSolution[i][j], 1.0);
        }
    }

    // Checking number of iterations within subiterations

    auto lvlOneMoisture = progressMoisture.getLevelOne();
    EXPECT_EQ(lvlOneMoisture, 650u);

    auto lvlTwoMoisture = progressMoisture.getLevelTwo();
    EXPECT_EQ(lvlTwoMoisture, 65000u);

    auto lvlThreeMoisture = progressMoisture.getLevelThree();
    EXPECT_EQ(lvlThreeMoisture, 0u);

    auto lvlOneThermal = progressThermal.getLevelOne();
    EXPECT_EQ(lvlOneThermal, 0u);

    auto lvlTwoThermal = progressThermal.getLevelTwo();
    EXPECT_EQ(lvlTwoThermal, 0u);

    auto lvlThreeThermal = progressThermal.getLevelThree();
    EXPECT_EQ(lvlThreeThermal, 0u);
}
