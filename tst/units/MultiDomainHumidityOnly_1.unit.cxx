#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::Nodes;

class MultiDomainHumidityOnly_1 : public testing::Test
{
protected:
    void SetUp() override
    {}

    void TearDown() override
    {}
};

//////////////////////////////////////////////////////////////////////////////////////////////////
// Test case with multi domain where only humidity calculations are performed while temperature is
// kept identical
//////////////////////////////////////////////////////////////////////////////////////////////////

TEST_F(MultiDomainHumidityOnly_1, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Simple two elements example with moisture and heat transfer.");

    const bool performThermalCalculations = false;
    const bool performMoistureCalculations = true;
    HygroThermFEM::MultiDomain multiDomain(performThermalCalculations, performMoistureCalculations);

    // Enter nodes. Arguments are: node number, x-coordinate, y-coordinate

    std::vector<double> gridXCoordinates{0, 0.05, 0.1};

    const double initialTemperature = 0.0;
    const double initialMoistureContent = 0.6;
    const double initialPressure = 101325;
    constexpr auto liquidPercent = 1.0;

    auto state = HygroThermFEM::State(
      initialTemperature, initialMoistureContent, initialPressure, liquidPercent);
    size_t nodeIndex = 0;
    for(auto val : gridXCoordinates)
    {
        ++nodeIndex;
        multiDomain.nodes().createNode(nodeIndex, val, 0.00, state);
        ++nodeIndex;
        multiDomain.nodes().createNode(nodeIndex, val, 0.05, state);
    }

    // Material Properties (Cottaer Sandstone)
    const auto & material = multiDomain.materials().createSolidMaterial({
        .name = "Cottaer Sandstone",
        .thermalConductivityDry = 1.8,
        .density = 2050.0,
        .porosity = 0.22,
        .heatCapacity = 850.0,
        .diffusionResistanceFactor = 15.0,
        .thermalConductivityMoistureDependent = {{0.0, 1.8}, {180, 1.8}},
        .moistureDependentMeasurementTemperature = 0,
        .thermalConductivityTemperatureDependent = {{0.0, 1.8}, {1, 1.8}},
        .temperatureDependentMeasurementHumidity = 0,
        .liquidTransportCurve = {{0, 0},
                                 {27, 1E-8},
                                 {45, 1.1E-8},
                                 {90, 2E-8},
                                 {126, 3.5E-8},
                                 {144, 5E-8},
                                 {162, 1E-7},
                                 {171, 2E-7},
                                 {180, 7E-7}},
        .sorptionCurve = {{0, 0},
                          {0.5, 5.3},
                          {0.65, 8.4},
                          {0.8, 12},
                          {0.93, 17},
                          {0.95, 25},
                          {0.99, 63},
                          {0.995, 83},
                          {0.999, 120},
                          {1, 180}}
    });

    /// Create elements
    for(size_t i = 1; i <= (multiDomain.nodes().maxIndex() - 2) / 2; ++i)
    {
        const auto node1 = 2u * i + 1u;
        const auto node2 = 2u * i + 2u;
        const auto node3 = 2u * i;
        const auto node4 = 2u * i - 1u;
        multiDomain.createElement(node1, node2, node3, node4, material.name());
    }

    /// Create Boundary Conditions
    constexpr auto hc = 1.0;
    constexpr auto airTemperature = 20.0;
    constexpr auto humidity = 0.0;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff{airTemperature, hc, humidity};

    multiDomain.createBC_FixedHc(1, 2, bcCoeff);

    constexpr auto dTime = 3600;
    constexpr auto nSteps = 10;

    auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
    auto humidities = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);
    std::vector<std::vector<double>> temperatureSolution;
    std::vector<std::vector<double>> waterContentSolution;

    for(auto i = 0; i < nSteps; ++i)
    {
        auto aSolution = multiDomain.transient(temperatures, humidities, dTime);
        temperatureSolution.push_back(aSolution.temperature);
        waterContentSolution.push_back(aSolution.waterContent);
        temperatures = aSolution.temperature;
        humidities = aSolution.humidity;
    }

    std::vector<std::vector<double>> correctWaterContentSolution{
      {7.055972, 7.055972, 7.366497, 7.366497, 7.366666, 7.366666},
      {6.748959, 6.748959, 7.366159, 7.366159, 7.366666, 7.366666},
      {6.444629, 6.444629, 7.365653, 7.365653, 7.366665, 7.366665},
      {6.143056, 6.143056, 7.364982, 7.364982, 7.366663, 7.366663},
      {5.844321, 5.844321, 7.364145, 7.364145, 7.366660, 7.366660},
      {5.548503, 5.548503, 7.363145, 7.363145, 7.366656, 7.366656},
      {5.277271, 5.277271, 7.361982, 7.361982, 7.366651, 7.366651},
      {5.128666, 5.128666, 7.360658, 7.360658, 7.366645, 7.366645},
      {4.981688, 4.981688, 7.359174, 7.359174, 7.366637, 7.366637},
      {4.836384, 4.836384, 7.357531, 7.357531, 7.366627, 7.366627}};

    EXPECT_EQ(waterContentSolution.size(), correctWaterContentSolution.size());

    for(auto i = 0u; i < correctWaterContentSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctWaterContentSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctWaterContentSolution[i][j], waterContentSolution[i][j], 1e-6);
        }
    }

    std::vector<std::vector<double>> correctTemperatureSolution{{0, 0, 0, 0, 0, 0},
                                                                {0, 0, 0, 0, 0, 0},
                                                                {0, 0, 0, 0, 0, 0},
                                                                {0, 0, 0, 0, 0, 0},
                                                                {0, 0, 0, 0, 0, 0},
                                                                {0, 0, 0, 0, 0, 0},
                                                                {0, 0, 0, 0, 0, 0},
                                                                {0, 0, 0, 0, 0, 0},
                                                                {0, 0, 0, 0, 0, 0},
                                                                {0, 0, 0, 0, 0, 0}};

    EXPECT_EQ(temperatureSolution.size(), correctTemperatureSolution.size());

    for(auto i = 0u; i < correctTemperatureSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctTemperatureSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctTemperatureSolution[i][j], temperatureSolution[i][j], 1e-6);
        }
    }
}
