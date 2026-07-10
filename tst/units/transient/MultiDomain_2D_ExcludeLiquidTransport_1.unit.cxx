#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "SlabCreator.hxx"
#include "TestHelpers.hxx"
#include "TestMaterials.hxx"

using HygroThermFEM::Nodes;
using HygroThermFEM::SimulationProperties;

class MultiDomain_2D_ExcludeLiquidTransport_1 : public testing::Test
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

TEST_F(MultiDomain_2D_ExcludeLiquidTransport_1, TestExample_1)
{
    constexpr auto excludeWaterLiquidTransportation{true};
    constexpr auto excludeHeatOfEvaporation{false};
    constexpr auto excludeCapillaryConduction{false};
    constexpr auto excludeVaporDiffusionConduction{false};
    constexpr auto thermalConductivityMoistureAndTemperatureDependent{false};

    SimulationProperties::Instance().setCalculationParameters(
      excludeWaterLiquidTransportation,
      excludeHeatOfEvaporation,
      excludeCapillaryConduction,
      excludeVaporDiffusionConduction,
      thermalConductivityMoistureAndTemperatureDependent);

    HygroThermFEM::MultiDomain multiDomain;

    const auto & material = multiDomain.materials().createSolidMaterial(TestHelper::CottaerSandstone());

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

    std::vector<std::vector<double>> correctWaterContentSolution{{60.3402861711, 60.3402861711, 62.9999422978, 62.9999422978, 60.3402861711, 60.3402861711},
 {57.6665849793, 57.6665849793, 62.999877525, 62.999877525, 57.6665849793, 57.6665849793},
 {54.9859133951, 54.9859133951, 62.9997527691, 62.9997527691, 54.9859133951, 54.9859133951},
 {52.3011756716, 52.3011756716, 62.9995608623, 62.9995608623, 52.3011756716, 52.3011756716},
 {49.6142625994, 49.6142625994, 62.9993017071, 62.9993017071, 49.6142625994, 49.6142625994},
 {46.9266144144, 46.9266144144, 62.9989761463, 62.9989761463, 46.9266144144, 46.9266144144},
 {44.2393660357, 44.2393660357, 62.998585018, 62.998585018, 44.2393660358, 44.2393660357},
 {41.5534161001, 41.5534161001, 62.9981290296, 62.9981290296, 41.5534161001, 41.5534161001},
 {38.8694826669, 38.8694826669, 62.9976087614, 62.9976087614, 38.8694826669, 38.8694826669},
 {36.1881349141, 36.1881349141, 62.9970246844, 62.9970246844, 36.1881349141, 36.1881349141},
 {33.509822764, 33.509822764, 62.9963771755, 62.9963771755, 33.509822764, 33.509822764},
 {30.8349011973, 30.8349011973, 62.9956665328, 62.9956665328, 30.8349011974, 30.8349011973},
 {28.1636497068, 28.1636497068, 62.9948929883, 62.9948929883, 28.1636497068, 28.1636497068},
 {25.4962878497, 25.4962878497, 62.9940567197, 62.9940567197, 25.4962878498, 25.4962878497},
 {22.8417960104, 22.8417960104, 62.993090614, 62.993090614, 22.8417960105, 22.8417960104},
 {20.2013287147, 20.2013287148, 62.9919815552, 62.9919815552, 20.2013287148, 20.2013287147},
 {17.5741551225, 17.5741551225, 62.9907298147, 62.9907298147, 17.5741551226, 17.5741551225},
 {15.0869737325, 15.0869737325, 62.9883648029, 62.9883648029, 15.0869737325, 15.0869737325},
 {12.7548875989, 12.7548875988, 62.9846871564, 62.9846871564, 12.7548875989, 12.7548875988},
 {10.6184111412, 10.6184111411, 62.9792868851, 62.9792868851, 10.6184111412, 10.6184111411},
 {8.68207589141, 8.68207589131, 62.9720020224, 62.9720020224, 8.68207589134, 8.68207589131},
 {6.9490124358, 6.9490124357, 62.9626756825, 62.9626756825, 6.94901243573, 6.94901243571},
 {5.40067280929, 5.40067280919, 62.9513126434, 62.9513126434, 5.40067280926, 5.4006728092},
 {4.16343502512, 4.16343502504, 62.9368121296, 62.9368121296, 4.1634350251, 4.16343502505}};

    TestHelper::dumpGolden("correctWaterContentSolution", waterContentSolution);
    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{{0.117508669586, 0.117508669586, 0.0808479736209, 0.0808479736209, 0.117508669586, 0.117508669586},
 {0.197289029523, 0.197289029523, 0.161007479652, 0.161007479652, 0.197289029523, 0.197289029523},
 {0.261539288222, 0.261539288222, 0.230241269327, 0.230241269327, 0.261539288222, 0.261539288222},
 {0.315722449476, 0.315722449476, 0.289134406496, 0.289134406496, 0.315722449476, 0.315722449476},
 {0.362260848034, 0.362260848034, 0.339542926439, 0.339542926439, 0.362260848034, 0.362260848034},
 {0.402787693615, 0.402787693615, 0.38317255143, 0.38317255143, 0.402787693615, 0.402787693615},
 {0.438561347675, 0.438561347675, 0.421422569674, 0.421422569674, 0.438561347675, 0.438561347675},
 {0.470539580501, 0.470539580501, 0.455375266867, 0.455375266867, 0.470539580501, 0.470539580501},
 {0.49951086615, 0.49951086615, 0.485919671016, 0.485919671016, 0.49951086615, 0.49951086615},
 {0.526106097089, 0.526106097089, 0.5137687345, 0.5137687345, 0.526106097089, 0.526106097089},
 {0.550827294737, 0.550827294737, 0.539490534236, 0.539490534236, 0.550827294737, 0.550827294737},
 {0.574073652639, 0.574073652639, 0.563537102894, 0.563537102894, 0.574073652639, 0.574073652639},
 {0.596162293822, 0.596162293822, 0.586267736164, 0.586267736164, 0.596162293822, 0.596162293822},
 {0.617444929831, 0.617444929831, 0.608071208089, 0.608071208089, 0.617444929831, 0.617444929831},
 {0.643239039081, 0.643239039081, 0.63272655914, 0.63272655914, 0.643239039081, 0.643239039081},
 {0.672913081269, 0.672913081269, 0.660949431251, 0.660949431251, 0.672913081269, 0.672913081269},
 {0.705341234228, 0.705341234228, 0.692186518097, 0.692186518097, 0.705341234228, 0.705341234228},
 {0.816324527945, 0.816324527945, 0.779351662201, 0.779351662201, 0.816324527944, 0.816324527945},
 {0.995411213987, 0.995411213988, 0.931407468105, 0.931407468105, 0.995411213987, 0.995411213987},
 {1.25751225665, 1.25751225665, 1.16191045519, 1.16191045519, 1.25751225665, 1.25751225665},
 {1.59395186305, 1.59395186306, 1.46926835164, 1.46926835164, 1.59395186305, 1.59395186306},
 {2.00017049528, 2.00017049528, 1.85017234026, 1.85017234026, 2.00017049528, 2.00017049528},
 {2.45957565225, 2.45957565226, 2.29207006668, 2.29207006668, 2.45957565225, 2.45957565226},
 {3.0553800585, 3.05538005851, 2.8516082316, 2.8516082316, 3.0553800585, 3.05538005851}};

    TestHelper::dumpGolden("correctTemperatureSolution", temperatureSolution);
    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-6);
}
