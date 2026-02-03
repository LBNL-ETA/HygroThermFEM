#include <memory>
#include <vector>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::NodePool;
using HygroThermFEM::MaterialPool;

class MultiDomain_2D_THERMM200mmSlab : public testing::Test
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

TEST_F(MultiDomain_2D_THERMM200mmSlab, TestExample_1)
{
    // Enter nodes. Arguments are: node number, x-coordinate, y-coordinate

    std::vector<double> gridXCoordinates{0, 0.05, 0.1};

    // const double initialTemperature = 0.0;
    // const double initialMoistureContent = 0.0;
    // const double initialPressure = 101325;
    // const auto liquidPercent = 1.0;

    // auto state = HygroThermFEM::State(
    //  initialTemperature, initialMoistureContent, initialPressure, liquidPercent);

    NodePool::Instance().createNode(1, 0.1, -0.049);
    NodePool::Instance().createNode(2, 0.1, -0.009);
    NodePool::Instance().createNode(3, 0.06, -0.049);
    NodePool::Instance().createNode(4, 0.1, 0.049);
    NodePool::Instance().createNode(5, 0.06, -0.009);
    NodePool::Instance().createNode(6, 0, -0.049);
    NodePool::Instance().createNode(7, 0.06, 0.049);
    NodePool::Instance().createNode(8, 0, -0.009);
    NodePool::Instance().createNode(9, -0.06, -0.049);
    NodePool::Instance().createNode(10, 0, 0.049);
    NodePool::Instance().createNode(11, -0.06, -0.009);
    NodePool::Instance().createNode(12, -0.1, -0.049);
    NodePool::Instance().createNode(13, -0.06, 0.049);
    NodePool::Instance().createNode(14, -0.1, -0.009);
    NodePool::Instance().createNode(15, -0.1, 0.049);

    // Material Properties
    constexpr double thermalConductivityDry{0.1};
    constexpr double density{400.0};
    constexpr double porosity{0.81};
    constexpr double specificHeatCapacityDry{850.0};
    constexpr double diffusionResistanceFactor{7.9};
    const std::vector<FenestrationCommon::point> thermalConductivityMoistureDependent = {
      {0.0, 0.1}, {380, 0.1}};
    constexpr double thermalConductivityMeasuredAtTemperature{0};
    const std::vector<FenestrationCommon::point> thermalConductivityTemperatureDependent = {
      {0.0, 0.1}, {1, 0.1}};
    constexpr double thermalConductivityMeasuredAtHumidity{0};
    const std::vector<FenestrationCommon::point> liquidTransportationCurve = {
      {0, 0}, {3.1, 2e-10}, {38, 5.4e-09}, {265, 1.1e-08}, {342, 2e-08}, {380, 1e-07}};

    const std::vector<FenestrationCommon::point> moistureStorageFunction = {{0, 0},
                                                                            {0.5, 3.1},
                                                                            {0.8, 8.4},
                                                                            {0.93, 18},
                                                                            {0.964, 116},
                                                                            {0.99, 266},
                                                                            {0.995, 269},
                                                                            {0.999, 277},
                                                                            {1, 380}};

    auto & material =
      MaterialPool::Instance().createSolidMaterial("Aerated Concrete (density: 400 kg/m)",
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
    domain.createElement(3, 1, 2, 5, material.name());
    domain.createElement(2, 4, 7, 5, material.name());
    domain.createElement(6, 3, 5, 8, material.name());
    domain.createElement(5, 7, 10, 8, material.name());
    domain.createElement(9, 6, 8, 11, material.name());
    domain.createElement(11, 14, 12, 9, material.name());
    domain.createElement(11, 8, 10, 13, material.name());
    domain.createElement(11, 13, 15, 14, material.name());

    /// Create Boundary Conditions
    const auto hc = 2.0;
    const auto airTemperature = 20.0;
    const auto humidity = 0.4;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff{airTemperature, hc, humidity};
    const std::vector<HygroThermFEM::FixedBCHCCoefficients> bcCoeffs{
      bcCoeff, bcCoeff, bcCoeff, bcCoeff, bcCoeff, bcCoeff, bcCoeff, bcCoeff, bcCoeff, bcCoeff,
      bcCoeff, bcCoeff, bcCoeff, bcCoeff, bcCoeff, bcCoeff, bcCoeff, bcCoeff, bcCoeff, bcCoeff};

    domain.createBC_FixedHc(15, 14, bcCoeffs);
    domain.createBC_FixedHc(14, 12, bcCoeffs);
    domain.createBC_FixedHc(1, 2, bcCoeffs);
    domain.createBC_FixedHc(2, 4, bcCoeffs);

    const auto dTime = 3600;
    const auto nSteps = 20;

    auto temperatures = NodePool::Instance().properties(HygroThermFEM::Variable::temperature);
    auto humidities = NodePool::Instance().properties(HygroThermFEM::Variable::humidity);
    std::vector<std::vector<double>> temperatureSolution;
    std::vector<std::vector<double>> humiditySolution;
    size_t timestepIndex{0u};

    for(auto i = 0; i < nSteps; ++i)
    {
        auto aSolution = domain.transient(temperatures, humidities, dTime, timestepIndex);
        temperatureSolution.push_back(aSolution.temperature);
        humiditySolution.push_back(aSolution.humidity);
        temperatures = aSolution.temperature;
        humidities = aSolution.humidity;
        ++timestepIndex;
    }

    // clang-format off
    std::vector<std::vector<double>> correctHumiditySolution{
        {0.212371,0.212371,0.001157,0.212371,0.001157,6e-06,0.001157,6e-06,0.001157,6e-06,0.001157,0.212371,0.001157,0.212371,0.212371},
        {0.317107,0.317107,0.003191,0.317107,0.003191,2.5e-05,0.003191,2.5e-05,0.003191,2.5e-05,0.003191,0.317107,0.003191,0.317107,0.317107},
        {0.373775,0.373775,0.00589,0.373775,0.00589,6.5e-05,0.00589,6.5e-05,0.00589,6.5e-05,0.00589,0.373775,0.00589,0.373775,0.373775},
        {0.40718,0.40718,0.009101,0.40718,0.009101,0.000134,0.009101,0.000134,0.009101,0.000134,0.009101,0.40718,0.009101,0.40718,0.40718},
        {0.427135,0.427135,0.012703,0.427135,0.012703,0.000242,0.012703,0.000242,0.012703,0.000242,0.012703,0.427135,0.012703,0.427135,0.427135},
        {0.438536,0.438536,0.016595,0.438536,0.016595,0.000396,0.016595,0.000396,0.016595,0.000396,0.016595,0.438536,0.016595,0.438536,0.438536},
        {0.444232,0.444232,0.020698,0.444232,0.020698,0.000603,0.020698,0.000603,0.020698,0.000603,0.020698,0.444232,0.020698,0.444232,0.444232},
        {0.446051,0.446051,0.024947,0.446051,0.024947,0.000868,0.024947,0.000868,0.024947,0.000868,0.024947,0.446051,0.024947,0.446051,0.446051},
        {0.445242,0.445242,0.02929,0.445242,0.02929,0.001195,0.02929,0.001195,0.02929,0.001195,0.02929,0.445242,0.02929,0.445242,0.445242},
        {0.442692,0.442692,0.033685,0.442692,0.033685,0.001588,0.033685,0.001588,0.033685,0.001588,0.033685,0.442692,0.033685,0.442692,0.442692},
        {0.439039,0.439039,0.038098,0.439039,0.038098,0.002049,0.038098,0.002049,0.038098,0.002049,0.038098,0.439039,0.038098,0.439039,0.439039},
        {0.434747,0.434747,0.042504,0.434747,0.042504,0.002577,0.042504,0.002577,0.042504,0.002577,0.042504,0.434747,0.042504,0.434747,0.434747},
        {0.430152,0.430152,0.046882,0.430152,0.046882,0.003175,0.046882,0.003175,0.046882,0.003175,0.046882,0.430152,0.046882,0.430152,0.430152},
        {0.425491,0.425491,0.051217,0.425491,0.051217,0.003841,0.051217,0.003841,0.051217,0.003841,0.051217,0.425491,0.051217,0.425491,0.425491},
        {0.420929,0.420929,0.055497,0.420929,0.055497,0.004574,0.055497,0.004574,0.055497,0.004574,0.055497,0.420929,0.055497,0.420929,0.420929},
        {0.416574,0.416574,0.059714,0.416574,0.059714,0.005374,0.059714,0.005374,0.059714,0.005374,0.059714,0.416574,0.059714,0.416574,0.416574},
        {0.412493,0.412493,0.063863,0.412493,0.063863,0.006239,0.063863,0.006239,0.063863,0.006239,0.063863,0.412493,0.063863,0.412493,0.412493},
        {0.408722,0.408722,0.067939,0.408722,0.067939,0.007166,0.067939,0.007166,0.067939,0.007166,0.067939,0.408722,0.067939,0.408722,0.408722},
        {0.405276,0.405276,0.07194,0.405276,0.07194,0.008155,0.07194,0.008155,0.07194,0.008155,0.07194,0.405276,0.07194,0.405276,0.405276},
        {0.402154,0.402154,0.075864,0.402154,0.075864,0.009202,0.075864,0.009202,0.075864,0.009202,0.075864,0.402154,0.075864,0.402154,0.402154}};
// clang-format on

    EXPECT_EQ(humiditySolution.size(), correctHumiditySolution.size());

    for(auto i = 0u; i < humiditySolution.size(); ++i)
    {
        for(auto j = 0u; j < humiditySolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctHumiditySolution[i][j], humiditySolution[i][j], 1e-6);
        }
    }

    // clang-format off
    std::vector<std::vector<double>> correctTemperatureSolution{
        {10.422203,10.422203,3.110194,10.422203,3.110194,1.151756,3.110194,1.151756,3.110194,1.151756,3.110194,10.422203,3.110194,10.422203,10.422203},
        {13.656363,13.656363,5.977403,13.656363,5.977403,2.938471,5.977403,2.938471,5.977403,2.938471,5.977403,13.656363,5.977403,13.656363,13.656363},
        {14.900768,14.900768,8.198622,14.900768,8.198622,4.885485,8.198622,4.885485,8.198622,4.885485,8.198622,14.900768,8.198622,14.900768,14.900768},
        {15.569701,15.569701,9.904939,15.569701,9.904939,6.742509,9.904939,6.742509,9.904939,6.742509,9.904939,15.569701,9.904939,15.569701,15.569701},
        {16.047302,16.047302,11.252232,16.047302,11.252232,8.409719,11.252232,8.409719,11.252232,8.409719,11.252232,16.047302,11.252232,16.047302,16.047302},
        {16.444317,16.444317,12.348744,16.444317,12.348744,9.864379,12.348744,9.864379,12.348744,9.864379,12.348744,16.444317,12.348744,16.444317,16.444317},
        {16.796594,16.796594,13.2634,16.796594,13.2634,11.117735,13.2634,11.117735,13.2634,11.117735,13.2634,16.796594,13.2634,16.796594,16.796594},
        {17.117655,17.117655,14.04056,17.117655,14.04056,12.193322,14.04056,12.193322,14.04056,12.193322,14.04056,17.117655,14.04056,17.117655,17.117655},
        {17.413474,17.413474,14.710001,17.413474,14.710001,13.117012,14.710001,13.117012,14.710001,13.117012,14.710001,17.413474,14.710001,17.413474,17.413474},
        {17.687101,17.687101,15.292674,17.687101,15.292674,13.912877,15.292674,13.912877,15.292674,13.912877,15.292674,17.687101,15.292674,17.687101,17.687101},
        {17.940324,17.940324,15.803945,17.940324,15.803945,14.601784,15.803945,14.601784,15.803945,14.601784,15.803945,17.940324,15.803945,17.940324,17.940324},
        {18.174354,18.174354,16.255455,18.174354,16.255455,15.2012,16.255455,15.2012,16.255455,15.2012,16.255455,18.174354,16.255455,18.174354,18.174354},
        {18.390137,18.390137,16.656239,18.390137,16.656239,15.725485,16.656239,15.725485,16.656239,15.725485,16.656239,18.390137,16.656239,18.390137,18.390137},
        {18.588551,18.588551,17.013455,18.588551,17.013455,16.186342,17.013455,16.186342,17.013455,16.186342,17.013455,18.588551,17.013455,18.588551,18.588551},
        {18.770428,18.770428,17.332861,18.770428,17.332861,16.593282,17.332861,16.593282,17.332861,16.593282,17.332861,18.770428,17.332861,18.770428,18.770428},
        {18.936626,18.936626,17.619158,18.936626,17.619158,16.954045,17.619158,16.954045,17.619158,16.954045,17.619158,18.936626,17.619158,18.936626,18.936626},
        {19.088035,19.088035,17.876245,19.088035,17.876245,17.274952,17.876245,17.274952,17.876245,17.274952,17.876245,19.088035,17.876245,19.088035,19.088035},
        {19.225571,19.225571,18.107398,19.225571,18.107398,17.561204,18.107398,17.561204,18.107398,17.561204,18.107398,19.225571,18.107398,19.225571,19.225571},
        {19.350169,19.350169,18.315412,19.350169,18.315412,17.817113,18.315412,17.817113,18.315412,17.817113,18.315412,19.350169,18.315412,19.350169,19.350169},
        {19.462765,19.462765,18.502699,19.462765,18.502699,18.046292,18.502699,18.046292,18.502699,18.046292,18.502699,19.462765,18.502699,19.462765,19.462765}};
    // clang-format on

    EXPECT_EQ(temperatureSolution.size(), correctTemperatureSolution.size());

    for(auto i = 0u; i < temperatureSolution.size(); ++i)
    {
        for(auto j = 0u; j < temperatureSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctTemperatureSolution[i][j], temperatureSolution[i][j], 1e-6);
        }
    }
}
