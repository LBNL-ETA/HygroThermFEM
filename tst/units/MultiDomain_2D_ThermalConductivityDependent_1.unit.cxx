#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "SlabCreator.hxx"
#include "TestMaterials.hxx"

using HygroThermFEM::Nodes;
using HygroThermFEM::SimulationProperties;

class MultiDomain_2D_ThermalConductivityDependent_1 : public testing::Test
{
protected:
    void SetUp() override
    {
        constexpr auto relaxationParameter{0.8};
        constexpr auto errorTolerance{1e-5};
        constexpr auto numberOfIterations{20u};
        SimulationProperties::Instance().setIterationParameters(
          relaxationParameter, errorTolerance, numberOfIterations);
    }

    void TearDown() override
    {
        SimulationProperties::Instance().reset();
    }
};

TEST_F(MultiDomain_2D_ThermalConductivityDependent_1, TestExample_1)
{
    constexpr auto excludeWaterLiquidTransportation{false};
    constexpr auto excludeHeatOfEvaporation{false};
    constexpr auto excludeCapillaryConduction{false};
    constexpr auto excludeVaporDiffusionConduction{false};
    constexpr auto thermalConductivityMoistureAndTemperatureDependent{true};

    SimulationProperties::Instance().setCalculationParameters(
      excludeWaterLiquidTransportation,
      excludeHeatOfEvaporation,
      excludeCapillaryConduction,
      excludeVaporDiffusionConduction,
      thermalConductivityMoistureAndTemperatureDependent);

    HygroThermFEM::MultiDomain multiDomain;

    auto params = TestHelper::CottaerSandstone();
    params.thermalConductivityMoistureDependent = {{0.0, 1.8}, {180, 2.5}};
    params.thermalConductivityTemperatureDependent = {{0.0, 1.8}, {1, 3.1}};
    const auto & material = multiDomain.materials().createSolidMaterial(params);

    constexpr HygroThermFEM::State state({
        .temperature = 0.0,
        .humidity = 0.99,
        .pressure = 101325.0,
        .liquidPercent = 1.0
    });

    TestHelper::SlabBuilder(multiDomain)
        .gridXCoordinates({0, 0.05, 0.1})
        .height(0.05)
        .material(material.name())
        .state(state)
        .startCorner(TestHelper::StartCorner::BottomRight)
        .direction(TestHelper::Direction::CounterClockwise)
        .build();

    /// Create Boundary Conditions
    constexpr auto hc = 5.0;
    constexpr auto airTemperature = 10.0;
    constexpr auto humidity = 0.0;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff{airTemperature, hc, humidity};

    multiDomain.createBC_FixedHc(1, 2, bcCoeff);
    multiDomain.createBC_FixedHc(5, 6, bcCoeff);

    constexpr auto dTime = 3600;
    constexpr auto nSteps = 24;

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
        {13.313629,13.313629,25.249588,25.249588,13.313629,13.313629},
        {12.614011,12.614011,25.243483,25.243483,12.614011,12.614011},
        {11.944162,11.944162,25.236386,25.236386,11.944162,11.944162},
        {11.503586,11.503586,25.228286,25.228286,11.503586,11.503586},
        {11.062114,11.062114,25.219166,25.219166,11.062114,11.062114},
        {10.620332,10.620332,25.209014,25.209014,10.620332,10.620332},
        {10.178957,10.178957,25.197814,25.197814,10.178957,10.178957},
        {9.742758,9.742758,25.185568,25.185568,9.742758,9.742758},
        {9.312013,9.312013,25.172273,25.172273,9.312013,9.312013},
        {8.886943,8.886943,25.157923,25.157923,8.886943,8.886943},
        {8.467731,8.467731,25.142515,25.142515,8.467731,8.467731},
        {8.102506,8.102506,25.126045,25.126045,8.102506,8.102506},
        {7.751963,7.751963,25.10851,25.10851,7.751963,7.751963},
        {7.40678,7.40678,25.089907,25.089907,7.40678,7.40678},
        {7.067029,7.067029,25.070232,25.070232,7.067029,7.067029},
        {6.732649,6.732649,25.04948,25.04948,6.732649,6.732649},
        {6.403487,6.403487,25.027648,25.027648,6.403487,6.403487},
        {6.079604,6.079604,25.00473,25.00473,6.079604,6.079604},
        {5.761057,5.761057,24.991884,24.991884,5.761057,5.761057},
        {5.447895,5.447895,24.981316,24.981316,5.447895,5.447895},
        {5.218018,5.218018,24.970287,24.970287,5.218018,5.218018},
        {5.062985,5.062985,24.958795,24.958795,5.062985,5.062985},
        {4.910774,4.910774,24.946838,24.946838,4.910774,4.910774},
        {4.761401,4.761401,24.934417,24.934417,4.761401,4.761401}};

    EXPECT_EQ(waterContentSolution.size(), correctWaterContentSolution.size());

    for(auto i = 0u; i < waterContentSolution.size(); ++i)
    {
        for(auto j = 0u; j < waterContentSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctWaterContentSolution[i][j], waterContentSolution[i][j], 1e-6);
        }
    }

    std::vector<std::vector<double>> correctTemperatureSolution{
        {0.346936,0.346936,0.269694,0.269694,0.346936,0.346936},
        {0.607102,0.607102,0.542994,0.542994,0.607102,0.607102},
        {0.835515,0.835515,0.786336,0.786336,0.835515,0.835515},
        {1.042046,1.042046,1.003074,1.003074,1.042046,1.042046},
        {1.23142,1.23142,1.199343,1.199343,1.23142,1.23142},
        {1.407397,1.407397,1.380142,1.380142,1.407397,1.407397},
        {1.573075,1.573075,1.549302,1.549302,1.573075,1.573075},
        {1.730663,1.730663,1.709507,1.709507,1.730663,1.730663},
        {1.881934,1.881934,1.862795,1.862795,1.881934,1.881934},
        {2.028277,2.028277,2.010729,2.010729,2.028277,2.028277},
        {2.170774,2.170774,2.154507,2.154507,2.170774,2.170774},
        {2.31026,2.31026,2.295025,2.295025,2.31026,2.31026},
        {2.447356,2.447356,2.433001,2.433001,2.447356,2.447356},
        {2.582566,2.582566,2.568961,2.568961,2.582566,2.582566},
        {2.716266,2.716266,2.703309,2.703309,2.716266,2.716266},
        {2.848748,2.848748,2.836359,2.836359,2.848748,2.848748},
        {2.980246,2.980246,2.96836,2.96836,2.980246,2.980246},
        {3.110933,3.110933,3.099498,3.099498,3.110933,3.110933},
        {3.240927,3.240927,3.229902,3.229902,3.240927,3.240927},
        {3.370314,3.370314,3.359663,3.359663,3.370314,3.370314},
        {3.499128,3.499128,3.488823,3.488823,3.499128,3.499128},
        {3.627386,3.627386,3.617405,3.617405,3.627386,3.627386},
        {3.75511,3.75511,3.745431,3.745431,3.75511,3.75511},
        {3.882303,3.882303,3.872908,3.872908,3.882303,3.882303}};

    EXPECT_EQ(temperatureSolution.size(), correctTemperatureSolution.size());

    for(auto i = 0u; i < temperatureSolution.size(); ++i)
    {
        for(auto j = 0u; j < temperatureSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctTemperatureSolution[i][j], temperatureSolution[i][j], 1e-6);
        }
    }
}
