#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

TEST(MultiDomain_2D_VariableTARPHc_MultiTimestepBC, TestExample_1)
{
    // Enter nodes. Arguments are: node number, x-coordinate, y-coordinate

    HygroThermFEM::MultiDomain multiDomain;

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
                        {1, 180}}});

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

    // Variable boundary conditions (temperature and humidity) over ten timesteps.
    const std::vector<HygroThermFEM::TARPCoefficients> bcCoeff{{20.0, 0.6},
                                                                       {20.0, 0.5},
                                                                       {20.0, 0.4},
                                                                       {20.0, 0.3},
                                                                       {20.0, 0.2},
                                                                       {18.0, 0.2},
                                                                       {16.0, 0.2},
                                                                       {14.0, 0.2},
                                                                       {12.0, 0.2},
                                                                       {10.0, 0.2}};

    const auto surfaceTilt{90.0};

    multiDomain.createBC_TARPHc(1, 2, bcCoeff, surfaceTilt);

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
      {2.962594, 2.962594, 0.003663, 0.003663, 8e-06, 8e-06},
      {4.409302, 4.409302, 0.009581, 0.009581, 3.4e-05, 3.4e-05},
      {5.131024, 5.131024, 0.016977, 0.016977, 8.4e-05, 8.4e-05},
      {5.452168, 5.452168, 0.025201, 0.025201, 0.000164, 0.000164},
      {5.254476, 5.254476, 0.033637, 0.033637, 0.000277, 0.000277},
      {5.059499, 5.059499, 0.042063, 0.042063, 0.000424, 0.000424},
      {4.811524, 4.811524, 0.050288, 0.050288, 0.000604, 0.000604},
      {4.527268, 4.527268, 0.058161, 0.058161, 0.000817, 0.000817},
      {4.220781, 4.220781, 0.065565, 0.065565, 0.001058, 0.001058},
      {3.903782, 3.903782, 0.072418, 0.072418, 0.001324, 0.001324}};

    EXPECT_EQ(waterContentSolution.size(), correctWaterContentSolution.size());

    for(auto i = 0u; i < correctWaterContentSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctWaterContentSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctWaterContentSolution[i][j], waterContentSolution[i][j], 1e-6);
        }
    }

    std::vector<std::vector<double>> correctTemperatureSolution{
      {3.696622, 3.696622, 1.920131, 1.920131, 1.437081, 1.437081},
      {5.359278, 5.359278, 3.642653, 3.642653, 3.087792, 3.087792},
      {6.376811, 6.376811, 4.989256, 4.989256, 4.510895, 4.510895},
      {7.075905, 7.075905, 6.009414, 6.009414, 5.632418, 5.632418},
      {7.563313, 7.563313, 6.766149, 6.766149, 6.48092, 6.48092},
      {7.909438, 7.909438, 7.321655, 7.321655, 7.110133, 7.110133},
      {8.125345, 8.125345, 7.710468, 7.710468, 7.559422, 7.559422},
      {8.218628, 8.218628, 7.953758, 7.953758, 7.854534, 7.854534},
      {8.200381, 8.200381, 8.06805, 8.06805, 8.014314, 8.014314},
      {8.080523, 8.080523, 8.066741, 8.066741, 8.05353, 8.05353}};

    EXPECT_EQ(temperatureSolution.size(), correctTemperatureSolution.size());

    for(auto i = 0u; i < correctTemperatureSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctTemperatureSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctTemperatureSolution[i][j], temperatureSolution[i][j], 1e-6);
        }        
    }
}
