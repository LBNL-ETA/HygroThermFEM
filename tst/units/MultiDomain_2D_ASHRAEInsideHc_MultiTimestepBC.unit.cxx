#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

TEST(MultiDomain_2D_ASHRAEInsideHc_MultiTimestepBC, TestExample_1)
{
    HygroThermFEM::MultiDomain multiDomain;

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

    // Variable boundary conditions (temperature, pressure and humidity) over ten timesteps.
    const std::vector<HygroThermFEM::ASHRAEInsideCoefficients> bcCoeff{{20.0, 0.6, 101325.0},
                                                                       {20.0, 0.5, 101325.0},
                                                                       {20.0, 0.4, 101325.0},
                                                                       {20.0, 0.3, 101325.0},
                                                                       {20.0, 0.2, 101325.0},
                                                                       {18.0, 0.2, 101325.0},
                                                                       {16.0, 0.2, 101325.0},
                                                                       {14.0, 0.2, 101325.0},
                                                                       {12.0, 0.2, 101325.0},
                                                                       {10.0, 0.2, 101325.0}};

    const auto surfaceTilt{90.0};    // degrees
    const auto surfaceHeight{1.0};   // meters

    multiDomain.createBC_ASHRAEInsideHc(1, 2, bcCoeff, surfaceHeight, surfaceTilt);

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
      {2.779608, 2.779608, 0.003408, 0.003408, 8e-06, 8e-06},
      {4.180509, 4.180509, 0.00893, 0.00893, 3.1e-05, 3.1e-05},
      {4.907596, 4.907596, 0.015857, 0.015857, 7.7e-05, 7.7e-05},
      {5.187632, 5.187632, 0.023601, 0.023601, 0.00015, 0.00015},
      {5.116003, 5.116003, 0.031604, 0.031604, 0.000253, 0.000253},
      {4.977132, 4.977132, 0.039663, 0.039663, 0.000387, 0.000387},
      {4.79477, 4.79477, 0.047614, 0.047614, 0.000553, 0.000553},
      {4.588012, 4.588012, 0.055336, 0.055336, 0.000748, 0.000748},
      {4.373407, 4.373407, 0.062748, 0.062748, 0.000972, 0.000972},
      {4.168545, 4.168545, 0.069817, 0.069817, 0.001221, 0.001221}};

    EXPECT_EQ(waterContentSolution.size(), correctWaterContentSolution.size());

    for(auto i = 0u; i < correctWaterContentSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctWaterContentSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctWaterContentSolution[i][j], waterContentSolution[i][j], 1e-6);
        }
    }

    std::vector<std::vector<double>> correctTemperatureSolution{
      {3.341972, 3.341972, 1.73596, 1.73596, 1.299243, 1.299243},
      {4.910432, 4.910432, 3.327175, 3.327175, 2.817003, 2.817003},
      {5.900677, 5.900677, 4.596291, 4.596291, 4.148667, 4.148667},
      {6.596693, 6.596693, 5.575781, 5.575781, 5.21675, 5.21675},
      {7.091264, 7.091264, 6.314986, 6.314986, 6.038688, 6.038688},
      {7.419941, 7.419941, 6.851804, 6.851804, 6.647232, 6.647232},
      {7.611637, 7.611637, 7.218784, 7.218784, 7.074981, 7.074981},
      {7.689664, 7.689664, 7.443674, 7.443674, 7.350903, 7.350903},
      {7.675293, 7.675293, 7.551012, 7.551012, 7.500651, 7.500651},
      {7.591928, 7.591928, 7.564889, 7.564889, 7.548708, 7.548708}};

    EXPECT_EQ(temperatureSolution.size(), correctTemperatureSolution.size());

    for(auto i = 0u; i < correctTemperatureSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctTemperatureSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctTemperatureSolution[i][j], temperatureSolution[i][j], 1e-6);
        }
    }
}
