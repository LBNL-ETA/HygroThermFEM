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
 {57.6665865749, 57.6665865749, 62.9998774972, 62.9998774972, 57.6665865749, 57.6665865749},
 {54.9859156652, 54.9859156652, 62.999752723, 62.999752723, 54.9859156652, 54.9859156652},
 {52.3011781241, 52.3011781241, 62.9995608032, 62.9995608032, 52.3011781241, 52.3011781241},
 {49.614264982, 49.614264982, 62.9993016381, 62.9993016381, 49.614264982, 49.614264982},
 {46.9266165712, 46.9266165712, 62.9989760695, 62.9989760695, 46.9266165712, 46.9266165712},
 {44.2393678596, 44.2393678596, 62.9985849353, 62.9985849353, 44.2393678596, 44.2393678596},
 {41.5534175157, 41.5534175157, 62.9981289429, 62.9981289429, 41.5534175157, 41.5534175157},
 {38.8694836215, 38.8694836215, 62.9976086723, 62.9976086723, 38.8694836215, 38.8694836215},
 {36.1881353718, 36.1881353718, 62.9970245945, 62.9970245945, 36.1881353719, 36.1881353718},
 {33.5098227007, 33.5098227007, 62.9963770865, 62.9963770865, 33.5098227007, 33.5098227007},
 {30.8349005969, 30.8349005969, 62.9956664461, 62.9956664461, 30.8349005969, 30.8349005969},
 {28.1636485585, 28.1636485585, 62.9948929055, 62.9948929055, 28.1636485585, 28.1636485585},
 {25.4962861461, 25.4962861461, 62.9940566423, 62.9940566423, 25.4962861461, 25.4962861461},
 {22.8417937565, 22.8417937565, 62.9930905433, 62.9930905433, 22.8417937565, 22.8417937565},
 {20.2013259203, 20.2013259203, 62.9919814924, 62.9919814924, 20.2013259203, 20.2013259203},
 {17.5741517373, 17.5741517373, 62.9907297617, 62.9907297617, 17.5741517373, 17.5741517373},
 {15.0869699482, 15.0869699482, 62.9883647601, 62.9883647601, 15.0869699482, 15.0869699482},
 {12.7548832084, 12.7548832084, 62.9846871245, 62.9846871245, 12.7548832084, 12.7548832084},
 {10.6184048968, 10.6184048968, 62.9792868776, 62.9792868776, 10.6184048968, 10.6184048968},
 {8.68206509007, 8.68206509005, 62.9720020702, 62.9720020702, 8.68206509017, 8.68206509014},
 {6.94899364624, 6.94899364622, 62.9626758339, 62.9626758339, 6.94899364634, 6.94899364631},
 {5.40064184264, 5.40064184265, 62.9513129725, 62.9513129725, 5.40064184278, 5.40064184272},
 {4.16339514504, 4.16339514505, 62.9368126836, 62.9368126836, 4.16339514515, 4.1633951451}};

    TestHelper::dumpGolden("correctWaterContentSolution", waterContentSolution);
    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{{0.117500136625, 0.117500136625, 0.0808591496619, 0.0808591496619, 0.117500136625, 0.117500136625},
 {0.197285440395, 0.197285440395, 0.161016806789, 0.161016806789, 0.197285440395, 0.197285440395},
 {0.261538332844, 0.261538332844, 0.230249412074, 0.230249412074, 0.261538332844, 0.261538332844},
 {0.3157228443, 0.3157228443, 0.28914178149, 0.28914178149, 0.3157228443, 0.3157228443},
 {0.362262078059, 0.362262078059, 0.339549564184, 0.339549564184, 0.362262078059, 0.362262078059},
 {0.402789498148, 0.402789498148, 0.383178391513, 0.383178391513, 0.402789498148, 0.402789498148},
 {0.438563558354, 0.438563558354, 0.421427550499, 0.421427550499, 0.438563558354, 0.438563558354},
 {0.470542076708, 0.470542076708, 0.45537934212, 0.45537934212, 0.470542076708, 0.470542076708},
 {0.49951355993, 0.49951355993, 0.485922808864, 0.485922808864, 0.49951355993, 0.49951355993},
 {0.526108925829, 0.526108925829, 0.513770913347, 0.513770913347, 0.526108925829, 0.526108925829},
 {0.550830215301, 0.550830215301, 0.539491739167, 0.539491739167, 0.550830215301, 0.550830215301},
 {0.574076636656, 0.574076636656, 0.563537322852, 0.563537322852, 0.574076636656, 0.574076636656},
 {0.596165323907, 0.596165323907, 0.586266961775, 0.586266961775, 0.596165323907, 0.596165323907},
 {0.617447998943, 0.617447998943, 0.608069429222, 0.608069429222, 0.617447998943, 0.617447998943},
 {0.643242086234, 0.643242086234, 0.632723948397, 0.632723948397, 0.643242086234, 0.643242086234},
 {0.672916436335, 0.672916436335, 0.660945811871, 0.660945811871, 0.672916436335, 0.672916436335},
 {0.70534514152, 0.70534514152, 0.692181675931, 0.692181675931, 0.70534514152, 0.70534514152},
 {0.816330280423, 0.816330280423, 0.779347340977, 0.779347340977, 0.816330280423, 0.816330280423},
 {0.995428822075, 0.995428822074, 0.931400425483, 0.931400425483, 0.995428822074, 0.995428822075},
 {1.25755615294, 1.25755615294, 1.16189649336, 1.16189649336, 1.25755615294, 1.25755615294},
 {1.59404138645, 1.59404138645, 1.46923968273, 1.46923968273, 1.59404138645, 1.59404138645},
 {2.00032767354, 2.00032767354, 1.85011674913, 1.85011674913, 2.00032767353, 2.00032767354},
 {2.45982013704, 2.45982013704, 2.29196983884, 2.29196983884, 2.45982013704, 2.45982013704},
 {3.05576378416, 3.05576378416, 2.85142749416, 2.85142749416, 3.05576378414, 3.05576378415}};

    TestHelper::dumpGolden("correctTemperatureSolution", temperatureSolution);
    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-6);
}
