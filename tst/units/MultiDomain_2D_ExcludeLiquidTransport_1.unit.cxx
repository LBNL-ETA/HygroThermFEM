#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::NodePool;
using HygroThermFEM::MaterialPool;
using HygroThermFEM::SimulationProperties;

class MultiDomain_2D_ExcludeLiquidTransport_1 : public testing::Test
{
protected:
    void SetUp() override
    {
        const auto relaxationParameter{0.8};
        const auto errorTolerance{1e-5};
        const auto numberOfIterations{20u};
        SimulationProperties::Instance().setIterationParameters(
          relaxationParameter, errorTolerance, numberOfIterations);
    }

    void TearDown() override
    {
        NodePool::Instance().clear();
        MaterialPool::Instance().clear();
        SimulationProperties::Instance().reset();
    }
};

TEST_F(MultiDomain_2D_ExcludeLiquidTransport_1, TestExample_1)
{
    const auto excludeWaterLiquidTransportation{true};
    const auto excludeHeatOfEvaporation{false};
    const auto excludeCapillaryConduction{false};
    const auto excludeVaporDiffusionConduction{false};
    const auto thermalConductivityMoistureAndTemperatureDependent{false};

    SimulationProperties::Instance().setCalculationParameters(
      excludeWaterLiquidTransportation,
      excludeHeatOfEvaporation,
      excludeCapillaryConduction,
      excludeVaporDiffusionConduction,
      thermalConductivityMoistureAndTemperatureDependent);

    std::vector<double> gridXCoordinates{0, 0.05, 0.1};

    const double initialTemperature = 0.0;
    const double initialMoistureContent = 0.99;
    const double initialPressure = 101325;
    const auto liquidPercent = 1.0;

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

    HygroThermFEM::MultiDomain domain;

    /// Create elements
    for(size_t i = 1; i <= (HygroThermFEM::maxNodeIndex() - 2) / 2; ++i)
    {
        const auto node1 = 2u * i + 1u;
        const auto node2 = 2u * i + 2u;
        const auto node3 = 2u * i;
        const auto node4 = 2u * i - 1u;
        createElement(domain, node1, node2, node3, node4, material.name());
    }

    /// Create Boundary Conditions
    const auto hc = 5.0;
    const auto airTemperature = 10.0;
    const auto humidity = 0.0;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff{airTemperature, hc, humidity};

    createBC_FixedHc(domain, 1, 2, bcCoeff);
    createBC_FixedHc(domain, 5, 6, bcCoeff);

    const auto dTime = 3600;
    const auto nSteps = 24;

    auto temperatures{properties(HygroThermFEM::Variable::temperature)};
    auto humidities{properties(HygroThermFEM::Variable::humidity)};
    std::vector<std::vector<double>> temperatureSolution;
    std::vector<std::vector<double>> waterContentSolution;

    HygroThermFEM::TransientSubstitutionSolver solver{domain};
    for(auto i = 0; i < nSteps; ++i)
    {
        auto aSolution = solver.transient(temperatures, humidities, dTime);
        temperatureSolution.push_back(aSolution.temperature);
        waterContentSolution.push_back(aSolution.waterContent);
        temperatures = aSolution.temperature;
        humidities = aSolution.humidity;
    }

    std::vector<std::vector<double>> correctWaterContentSolution{
      {62.286722, 62.286722, 63.000012, 63.000012, 62.286722, 62.286722},
      {61.552789, 61.552789, 63.000003, 63.000003, 61.552789, 61.552789},
      {60.797261, 60.797261, 62.999993, 62.999993, 60.797261, 60.797261},
      {60.018854, 60.018854, 62.999978, 62.999978, 60.018854, 60.018854},
      {59.215961, 59.215961, 62.999957, 62.999957, 59.215961, 59.215961},
      {58.386612, 58.386612, 62.999930, 62.999930, 58.386612, 58.386612},
      {57.528416, 57.528416, 62.999896, 62.999896, 57.528416, 57.528416},
      {56.638464, 56.638464, 62.999856, 62.999856, 56.638464, 56.638464},
      {55.713204, 55.713204, 62.999808, 62.999808, 55.713204, 55.713204},
      {54.748254, 54.748254, 62.999753, 62.999753, 54.748254, 54.748254},
      {53.738140, 53.738140, 62.999691, 62.999691, 53.738140, 53.738140},
      {52.675896, 52.675896, 62.999620, 62.999620, 52.675896, 52.675896},
      {51.552443, 51.552443, 62.999540, 62.999540, 51.552443, 51.552443},
      {50.355594, 50.355594, 62.999452, 62.999452, 50.355594, 50.355594},
      {49.068328, 49.068328, 62.999352, 62.999352, 49.068328, 49.068328},
      {47.665556, 47.665556, 62.999242, 62.999242, 47.665556, 47.665556},
      {46.107403, 46.107403, 62.999119, 62.999119, 46.107403, 46.107403},
      {44.322661, 44.322661, 62.998980, 62.998980, 44.322661, 44.322661},
      {42.455999, 42.455999, 62.998827, 62.998827, 42.455999, 42.455999},
      {40.549330, 40.549330, 62.998659, 62.998659, 40.549330, 40.549330},
      {38.599831, 38.599831, 62.998477, 62.998477, 38.599831, 38.599831},
      {36.604325, 36.604325, 62.998278, 62.998278, 36.604325, 36.604325},
      {34.559222, 34.559222, 62.998063, 62.998063, 34.559222, 34.559222},
      {32.460437, 32.460437, 62.997832, 62.997832, 32.460437, 32.460437},
    };

    EXPECT_EQ(waterContentSolution.size(), correctWaterContentSolution.size());

    for(auto i = 0u; i < waterContentSolution.size(); ++i)
    {
        for(auto j = 0u; j < waterContentSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctWaterContentSolution[i][j], waterContentSolution[i][j], 1e-6);
        }
    }

    std::vector<std::vector<double>> correctTemperatureSolution{
      {0.103036, 0.103036, 0.070870, 0.070870, 0.103036, 0.103036},
      {0.173784, 0.173784, 0.141661, 0.141661, 0.173784, 0.173784},
      {0.229594, 0.229594, 0.202149, 0.202149, 0.229594, 0.229594},
      {0.275310, 0.275310, 0.252476, 0.252476, 0.275310, 0.275310},
      {0.313224, 0.313224, 0.294262, 0.294262, 0.313224, 0.313224},
      {0.344908, 0.344908, 0.329097, 0.329097, 0.344908, 0.344908},
      {0.371589, 0.371589, 0.358321, 0.358321, 0.371589, 0.371589},
      {0.394260, 0.394260, 0.383034, 0.383034, 0.394260, 0.394260},
      {0.413728, 0.413728, 0.404138, 0.404138, 0.413728, 0.413728},
      {0.430655, 0.430655, 0.422367, 0.422367, 0.430655, 0.430655},
      {0.445584, 0.445584, 0.438326, 0.438326, 0.445584, 0.445584},
      {0.458965, 0.458965, 0.452514, 0.452514, 0.458965, 0.458965},
      {0.471176, 0.471176, 0.465345, 0.465345, 0.471176, 0.471176},
      {0.482539, 0.482539, 0.477171, 0.477171, 0.482539, 0.482539},
      {0.493339, 0.493339, 0.488298, 0.488298, 0.493339, 0.493339},
      {0.503841, 0.503841, 0.499006, 0.499006, 0.503841, 0.503841},
      {0.514317, 0.514317, 0.509568, 0.509568, 0.514317, 0.514317},
      {0.525088, 0.525088, 0.520294, 0.520294, 0.525088, 0.525088},
      {0.536156, 0.536156, 0.531269, 0.531269, 0.536156, 0.536156},
      {0.547500, 0.547500, 0.542514, 0.542514, 0.547500, 0.547500},
      {0.559129, 0.559129, 0.554038, 0.554038, 0.559129, 0.559129},
      {0.571052, 0.571052, 0.565856, 0.565856, 0.571052, 0.571052},
      {0.583288, 0.583288, 0.577982, 0.577982, 0.583288, 0.583288},
      {0.595854, 0.595854, 0.590434, 0.590434, 0.595854, 0.595854},
    };

    EXPECT_EQ(temperatureSolution.size(), correctTemperatureSolution.size());

    for(auto i = 0u; i < temperatureSolution.size(); ++i)
    {
        for(auto j = 0u; j < temperatureSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctTemperatureSolution[i][j], temperatureSolution[i][j], 1e-6);
        }
    }
}
