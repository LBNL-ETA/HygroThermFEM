#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "PrintResults.hxx"

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

class MultiDomain_HighHumidity : public testing::Test
{
public:
    HygroThermFEM::MultiDomain domain;

    const double dTime{360};
    const size_t nSteps{5u};

    const double domainTemperature{0.0};
    const double domainHumidity{0.999};
    const double domainPressure{101325.0};
    const double liquidPercent{1.0};

    // Create Boundary Conditions
    const double hc{10.0};
    const double airTemperature{20.0};
    const double airHumidity{1.0};

protected:
    void SetUp() override
    {
        std::vector<double> gridXCoordinates{0.15, 0.05, 0.00};

        HygroThermFEM::State state(
          domainTemperature, domainHumidity, domainPressure, liquidPercent);
        size_t nodeIndex = 0;
        for(auto val : gridXCoordinates)
        {
            ++nodeIndex;
            NodePool::Instance().createNode(nodeIndex, val, 0.05, state);
            ++nodeIndex;
            NodePool::Instance().createNode(nodeIndex, val, 0.00, state);
        }

        // Material Properties (Cottaer Sandstone)
        const double thermalConductivityDry{1.8};
        const double density{2050.0};
        const double porosity{0.22};
        const double specificHeatCapacityDry{850.0};
        const double diffusionResistanceFactor{15.0};
        const std::vector<FenestrationCommon::point> thermalConductivityMoistureDependent = {
          {0.0, 1.8}, {180, 1.8}};
        const double thermalConductivityMeasuredAtTemperature{0};
        const std::vector<FenestrationCommon::point> thermalConductivityTemperatureDependent = {
          {0.0, 1.8}, {1, 1.8}};
        const double thermalConductivityMeasuredAtHumidity{0};
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

        /// Create elements
        for(size_t i = 1; i <= (HygroThermFEM::maxNodeIndex() - 2) / 2; ++i)
        {
            const auto node1 = 2u * i + 1u;
            const auto node2 = 2u * i + 2u;
            const auto node3 = 2u * i;
            const auto node4 = 2u * i - 1u;
            createElement(domain, node2, node3, node4, node1, material.name());
        }

        const HygroThermFEM::FixedBCHCCoefficients bcCoeff{airTemperature, hc, airHumidity};

        createBC_FixedHc(domain, 5, 6, bcCoeff);
    }

    void TearDown() override
    {
        NodePool::Instance().clear();
        MaterialPool::Instance().clear();
    }
};

TEST_F(MultiDomain_HighHumidity, Substitution)
{
    SCOPED_TRACE("Begin Test: Simple two elements example with moisture transfer.");

    ObserveSimulationProgrees progressThermal;
    subscribeThermal(domain, &progressThermal);

    ObserveSimulationProgrees progressMoisture;
    subscribeMoisture(domain, &progressMoisture);

    auto temperatures{properties(HygroThermFEM::Variable::temperature)};
    auto humidities{properties(HygroThermFEM::Variable::humidity)};
    std::vector<std::vector<double>> temperatureSolution;
    std::vector<double> temperatureError;
    std::vector<std::vector<double>> waterContentSolution;
    std::vector<double> humidityError;
    size_t timestepIndex{0};

    HygroThermFEM::TransientSubstitutionSolver solver{domain};
    for(auto i = 0u; i < nSteps; ++i)
    {
        auto aSolution = solver.transient(temperatures, humidities, dTime, timestepIndex);
        temperatureSolution.push_back(aSolution.temperature);
        temperatureError.push_back(aSolution.temperatureError);
        waterContentSolution.push_back(aSolution.waterContent);
        humidityError.push_back(aSolution.humidityError);
        temperatures = aSolution.temperature;
        humidities = aSolution.humidity;
        ++timestepIndex;
    }

    const std::vector<double> correctHumidityError{
      7.891482e-07, 2.315632e-07, 1.651598e-07, 1.257494e-06, 8.297355e-06};
    const std::vector<std::vector<double>> correctWaterContentSolution{
      {120.131097, 120.131097, 120.277551, 120.277551, 120.486078, 120.486078},
      {120.373681, 120.373681, 120.534680, 120.534680, 120.736459, 120.736459},
      {120.619012, 120.619012, 120.774788, 120.774788, 120.967641, 120.967641},
      {120.863364, 120.863364, 120.899965, 120.899965, 121.487909, 121.487909},
      {120.915509, 120.915509, 121.971323, 121.971323, 119.924717, 119.924717}};

    EXPECT_EQ(waterContentSolution.size(), correctWaterContentSolution.size());

    for(auto i = 0u; i < waterContentSolution.size(); ++i)
    {
        EXPECT_NEAR(correctHumidityError[i], humidityError[i], 1e-10);
        for(auto j = 0u; j < waterContentSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctWaterContentSolution[i][j], waterContentSolution[i][j], 1e-6);
        }
    }

    const std::vector<double> correctTemperatureError{
      2.067586e-04, 3.417621e-04, 2.755778e-04, 9.723013e-05, 6.292538e-04};
    const std::vector<std::vector<double>> correctTemperatureSolution{
      {0.003947, 0.003947, 0.093358, 0.093358, 1.724834, 1.724834},
      {0.014523, 0.014523, 0.254303, 0.254303, 3.112191, 3.112191},
      {0.033466, 0.033466, 0.463318, 0.463318, 4.237048, 4.237048},
      {0.061819, 0.061819, 0.705711, 0.705711, 5.157082, 5.157082},
      {0.100096, 0.100096, 0.970019, 0.970019, 5.921223, 5.921223}};

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
    EXPECT_EQ(lvlOneMoisture, 5u);

    auto lvlTwoMoisture = progressMoisture.getLevelTwo();
    EXPECT_EQ(lvlTwoMoisture, 500u);

    auto lvlThreeMoisture = progressMoisture.getLevelThree();
    EXPECT_EQ(lvlThreeMoisture, 0u);

    auto lvlOneThermal = progressThermal.getLevelOne();
    EXPECT_EQ(lvlOneThermal, 0u);

    auto lvlTwoThermal = progressThermal.getLevelTwo();
    EXPECT_EQ(lvlTwoThermal, 0u);

    auto lvlThreeThermal = progressThermal.getLevelThree();
    EXPECT_EQ(lvlThreeThermal, 0u);
}

TEST_F(MultiDomain_HighHumidity, Sundials)
{
    SCOPED_TRACE("Begin Test: Simple two elements example with moisture transfer.");

    auto temperatures{properties(HygroThermFEM::Variable::temperature)};
    auto humidities{properties(HygroThermFEM::Variable::humidity)};
    std::vector<std::vector<double>> temperatureSolution;
    std::vector<double> temperatureError;
    std::vector<std::vector<double>> waterContentSolution;
    std::vector<double> humidityError;
    size_t timestepIndex{0};

    Sundials::SolverIDA solver{domain};
    for(auto i = 0u; i < nSteps; ++i)
    {
        auto aSolution = solver.transient(temperatures, humidities, dTime, i);
        temperatureSolution.push_back(aSolution.temperature);
        temperatureError.push_back(aSolution.temperatureError);
        waterContentSolution.push_back(aSolution.waterContent);
        humidityError.push_back(aSolution.humidityError);
        temperatures = aSolution.temperature;
        humidities = aSolution.humidity;
        ++timestepIndex;
    }

    const std::vector<double> correctHumidityError{
      0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00};
    const std::vector<std::vector<double>> correctWaterContentSolution{
      {120.000002, 120.000002, 120.002660, 120.002660, 121.695731, 722.295140},
      {120.232754, 121.766152, 1214.566794, 669.942953, 1708.520796, 16.266035},
      {2343.596285, 2325.722773, 46.933135, 49.542321, 57.681948, 32689.839057},
      {1110.148110, 1206.699228, 11547.899929, 6388.912307, 12008.837066, 19765.411391},
      {1510.028538, 2168.949978, 10775.578976, 7424.187594, 11433.807370, 16635.962366}};

    EXPECT_EQ(waterContentSolution.size(), correctWaterContentSolution.size());

    for(auto i = 0u; i < waterContentSolution.size(); ++i)
    {
        EXPECT_NEAR(correctHumidityError[i], humidityError[i], 1e-10);
        for(auto j = 0u; j < waterContentSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctWaterContentSolution[i][j], waterContentSolution[i][j], 1e-6);
        }
    }

    const std::vector<double> correctTemperatureError{
      0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00};
    const std::vector<std::vector<double>> correctTemperatureSolution{
      {0.000682, 0.000295, 0.022362, 0.032322, 1.545948, 0.340728},
      {0.002419, 0.001994, 0.037004, 0.043627, 1.756493, 1.453173},
      {0.002250, 0.001841, 0.053068, 0.057039, 1.777646, 1.462069},
      {0.002835, 0.002813, 0.055405, 0.063991, 1.783232, 1.457015},
      {0.003397, 0.003341, 0.057253, 0.067956, 1.788020, 1.461584}};

    EXPECT_EQ(temperatureSolution.size(), correctTemperatureSolution.size());

    for(auto i = 0u; i < correctTemperatureSolution.size(); ++i)
    {
        EXPECT_NEAR(correctTemperatureError[i], temperatureError[i], 1e-10);
        for(auto j = 0u; j < correctTemperatureSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctTemperatureSolution[i][j], temperatureSolution[i][j], 1e-6);
        }
    }
}
