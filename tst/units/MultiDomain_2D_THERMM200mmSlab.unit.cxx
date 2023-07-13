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
    const double thermalConductivityDry{0.1};
    const double density{400.0};
    const double porosity{0.81};
    const double specificHeatCapacityDry{850.0};
    const double diffusionResistanceFactor{7.9};
    const std::vector<FenestrationCommon::point> thermalConductivityMoistureDependent = {
      {0.0, 0.1}, {380, 0.1}};
    const double thermalConductivityMeasuredAtTemperature{0};
    const std::vector<FenestrationCommon::point> thermalConductivityTemperatureDependent = {
      {0.0, 0.1}, {1, 0.1}};
    const double thermalConductivityMeasuredAtHumidity{0};
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
    createElement(domain, 3, 1, 2, 5, material.name());
    createElement(domain, 2, 4, 7, 5, material.name());
    createElement(domain, 6, 3, 5, 8, material.name());
    createElement(domain, 5, 7, 10, 8, material.name());
    createElement(domain, 9, 6, 8, 11, material.name());
    createElement(domain, 11, 14, 12, 9, material.name());
    createElement(domain, 11, 8, 10, 13, material.name());
    createElement(domain, 11, 13, 15, 14, material.name());

    /// Create Boundary Conditions
    const auto hc = 2.0;
    const auto airTemperature = 20.0;
    const auto humidity = 0.4;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff{airTemperature, hc, humidity};
    const std::vector<HygroThermFEM::FixedBCHCCoefficients> bcCoeffs{
      bcCoeff, bcCoeff, bcCoeff, bcCoeff, bcCoeff, bcCoeff, bcCoeff, bcCoeff, bcCoeff, bcCoeff,
      bcCoeff, bcCoeff, bcCoeff, bcCoeff, bcCoeff, bcCoeff, bcCoeff, bcCoeff, bcCoeff, bcCoeff};

    createBC_FixedHc(domain, 15, 14, bcCoeffs);
    createBC_FixedHc(domain, 14, 12, bcCoeffs);
    createBC_FixedHc(domain, 1, 2, bcCoeffs);
    createBC_FixedHc(domain, 2, 4, bcCoeffs);

    const auto dTime = 3600;
    const auto nSteps = 20;

    auto temperatures{properties(HygroThermFEM::Variable::temperature)};
    auto humidities{properties(HygroThermFEM::Variable::humidity)};
    std::vector<std::vector<double>> temperatureSolution;
    std::vector<std::vector<double>> humiditySolution;
    size_t timestepIndex{0u};

    for(auto i = 0; i < nSteps; ++i)
    {
        HygroThermFEM::TransientSubstitutionSolver solver;
        auto aSolution = solver.transient(domain, temperatures, humidities, dTime, timestepIndex);
        temperatureSolution.push_back(aSolution.temperature);
        humiditySolution.push_back(aSolution.humidity);
        temperatures = aSolution.temperature;
        humidities = aSolution.humidity;
        ++timestepIndex;
    }

    // clang-format off
    std::vector<std::vector<double>> correctHumiditySolution{
        {0.212372,0.212372,0.001157,0.212372,0.001157,0.000006,0.001157,0.000006,0.001157,0.000006,0.001157,0.212372,0.001157,0.212372,0.212372},
        {0.317107,0.317107,0.003191,0.317107,0.003191,0.000025,0.003191,0.000025,0.003191,0.000025,0.003191,0.317107,0.003191,0.317107,0.317107},
        {0.373775,0.373775,0.005890,0.373775,0.005890,0.000065,0.005890,0.000065,0.005890,0.000065,0.005890,0.373775,0.005890,0.373775,0.373775},
        {0.407179,0.407179,0.009102,0.407179,0.009102,0.000135,0.009102,0.000135,0.009102,0.000135,0.009102,0.407179,0.009102,0.407179,0.407179},
        {0.427134,0.427134,0.012703,0.427134,0.012703,0.000243,0.012703,0.000243,0.012703,0.000243,0.012703,0.427134,0.012703,0.427134,0.427134},
        {0.438535,0.438535,0.016595,0.438535,0.016595,0.000397,0.016595,0.000397,0.016595,0.000397,0.016595,0.438535,0.016595,0.438535,0.438535},
        {0.444230,0.444230,0.020698,0.444230,0.020698,0.000604,0.020698,0.000604,0.020698,0.000604,0.020698,0.444230,0.020698,0.444230,0.444230},
        {0.446049,0.446049,0.024948,0.446049,0.024948,0.000869,0.024948,0.000869,0.024948,0.000869,0.024948,0.446049,0.024948,0.446049,0.446049},
        {0.445241,0.445241,0.029291,0.445241,0.029291,0.001196,0.029291,0.001196,0.029291,0.001196,0.029291,0.445241,0.029291,0.445241,0.445241},
        {0.442690,0.442690,0.033685,0.442690,0.033685,0.001589,0.033685,0.001589,0.033685,0.001589,0.033685,0.442690,0.033685,0.442690,0.442690},
        {0.439037,0.439037,0.038099,0.439037,0.038099,0.002049,0.038099,0.002049,0.038099,0.002049,0.038099,0.439037,0.038099,0.439037,0.439037},
        {0.434746,0.434746,0.042505,0.434746,0.042505,0.002578,0.042505,0.002578,0.042505,0.002578,0.042505,0.434746,0.042505,0.434746,0.434746},
        {0.430153,0.430153,0.046883,0.430153,0.046883,0.003175,0.046883,0.003175,0.046883,0.003175,0.046883,0.430153,0.046883,0.430153,0.430153},
        {0.425494,0.425494,0.051217,0.425494,0.051217,0.003841,0.051217,0.003841,0.051217,0.003841,0.051217,0.425494,0.051217,0.425494,0.425494},
        {0.420933,0.420933,0.055498,0.420933,0.055498,0.004575,0.055498,0.004575,0.055498,0.004575,0.055498,0.420933,0.055498,0.420933,0.420933},
        {0.416578,0.416578,0.059715,0.416578,0.059715,0.005375,0.059715,0.005375,0.059715,0.005375,0.059715,0.416578,0.059715,0.416578,0.416578},
        {0.412497,0.412497,0.063863,0.412497,0.063863,0.006239,0.063863,0.006239,0.063863,0.006239,0.063863,0.412497,0.063863,0.412497,0.412497},
        {0.408727,0.408727,0.067939,0.408727,0.067939,0.007167,0.067939,0.007167,0.067939,0.007167,0.067939,0.408727,0.067939,0.408727,0.408727},
        {0.405280,0.405280,0.071940,0.405280,0.071940,0.008155,0.071940,0.008155,0.071940,0.008155,0.071940,0.405280,0.071940,0.405280,0.405280},
        {0.402158,0.402158,0.075865,0.402158,0.075865,0.009203,0.075865,0.009203,0.075865,0.009203,0.075865,0.402158,0.075865,0.402158,0.402158},
    };
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
        {13.656355,13.656355,5.977401,13.656355,5.977401,2.938471,5.977401,2.938471,5.977401,2.938471,5.977401,13.656355,5.977401,13.656355,13.656355},
        {14.900757,14.900757,8.198617,14.900757,8.198617,4.885483,8.198617,4.885483,8.198617,4.885483,8.198617,14.900757,8.198617,14.900757,14.900757},
        {15.569689,15.569689,9.904932,15.569689,9.904932,6.742505,9.904932,6.742505,9.904932,6.742505,9.904932,15.569689,9.904932,15.569689,15.569689},
        {16.047290,16.047290,11.252223,16.047290,11.252223,8.409713,11.252223,8.409713,11.252223,8.409713,11.252223,16.047290,11.252223,16.047290,16.047290},
        {16.444305,16.444305,12.348734,16.444305,12.348734,9.864372,12.348734,9.864372,12.348734,9.864372,12.348734,16.444305,12.348734,16.444305,16.444305},
        {16.796582,16.796582,13.263390,16.796582,13.263390,11.117727,13.263390,11.117727,13.263390,11.117727,13.263390,16.796582,13.263390,16.796582,16.796582},
        {17.117644,17.117644,14.040550,17.117644,14.040550,12.193313,14.040550,12.193313,14.040550,12.193313,14.040550,17.117644,14.040550,17.117644,17.117644},
        {17.413465,17.413465,14.709992,17.413465,14.709992,13.117003,14.709992,13.117003,14.709992,13.117003,14.709992,17.413465,14.709992,17.413465,17.413465},
        {17.687095,17.687095,15.292665,17.687095,15.292665,13.912868,15.292665,13.912868,15.292665,13.912868,15.292665,17.687095,15.292665,17.687095,17.687095},
        {17.940320,17.940320,15.803938,17.940320,15.803938,14.601776,15.803938,14.601776,15.803938,14.601776,15.803938,17.940320,15.803938,17.940320,17.940320},
        {18.174352,18.174352,16.255449,18.174352,16.255449,15.201193,16.255449,15.201193,16.255449,15.201193,16.255449,18.174352,16.255449,18.174352,18.174352},
        {18.390197,18.390197,16.656255,18.390197,16.656255,15.725486,16.656255,15.725486,16.656255,15.725486,16.656255,18.390197,16.656255,18.390197,18.390197},
        {18.588611,18.588611,17.013484,18.588611,17.013484,16.186353,17.013484,16.186353,17.013484,16.186353,17.013484,18.588611,17.013484,18.588611,18.588611},
        {18.770477,18.770477,17.332894,18.770477,17.332894,16.593302,17.332894,16.593302,17.332894,16.593302,17.332894,18.770477,17.332894,18.770477,18.770477},
        {18.936662,18.936662,17.619190,18.936662,17.619190,16.954069,17.619190,16.954069,17.619190,16.954069,17.619190,18.936662,17.619190,18.936662,18.936662},
        {19.088060,19.088060,17.876274,19.088060,17.876274,17.274978,17.876274,17.274978,17.876274,17.274978,17.876274,19.088060,17.876274,19.088060,19.088060},
        {19.225587,19.225587,18.107423,19.225587,18.107423,17.561230,18.107423,17.561230,18.107423,17.561230,18.107423,19.225587,18.107423,19.225587,19.225587},
        {19.350178,19.350178,18.315432,19.350178,18.315432,17.817137,18.315432,17.817137,18.315432,17.817137,18.315432,19.350178,18.315432,19.350178,19.350178},
        {19.462769,19.462769,18.502714,19.462769,18.502714,18.046313,18.502714,18.046313,18.502714,18.046313,18.502714,19.462769,18.502714,19.462769,19.462769},
    };
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
