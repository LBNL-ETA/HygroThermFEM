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

    std::vector<std::vector<double>> correctWaterContentSolution{{60.3402861712, 60.3402861712, 62.9999455855, 62.9999455855, 60.3402861712, 60.3402861712},
 {57.6665850129, 57.6665850129, 62.999884098, 62.999884098, 57.6665850129, 57.6665850129},
 {54.9859120182, 54.9859120182, 62.999762512, 62.999762512, 54.9859120182, 54.9859120182},
 {52.3011767023, 52.3011767023, 62.9995736174, 62.9995736174, 52.3011767023, 52.3011767023},
 {49.6142781339, 49.6142781339, 62.9993176832, 62.9993176832, 49.6142781339, 49.6142781339},
 {46.9266672632, 46.9266672632, 62.9989921372, 62.9989921372, 46.9266672632, 46.9266672632},
 {44.2394916532, 44.2394916527, 62.9986010566, 62.9986010566, 44.2394916523, 44.2394916527},
 {41.5536638668, 41.5536638668, 62.9981451671, 62.9981451671, 41.5536638659, 41.5536638664},
 {38.8699093438, 38.8699093433, 62.9976250641, 62.9976250641, 38.8699093433, 38.8699093429},
 {36.188803991, 36.1888039901, 62.9970412308, 62.9970412303, 36.188803991, 36.1888039896},
 {33.5108043541, 33.5108043541, 62.9963940578, 62.9963940565, 33.5108043546, 33.5108043537},
 {30.8362718696, 30.83627187, 62.9956838562, 62.9956838549, 30.83627187, 30.83627187},
 {28.1654923219, 28.1654923219, 62.9949108728, 62.9949108711, 28.1654923219, 28.1654923219},
 {25.4986914169, 25.4986914169, 62.9940752994, 62.9940752977, 25.4986914169, 25.4986914169},
 {22.8448549521, 22.8448549516, 62.9931101651, 62.9931101633, 22.8448549518, 22.8448549516},
 {20.2051892362, 20.2051892346, 62.9920023357, 62.992002334, 20.2051892349, 20.2051892347},
 {17.57898521, 17.5789852085, 62.9907521582, 62.9907521585, 17.5789852088, 17.5789852085},
 {15.0925357023, 15.0925356992, 62.9883924161, 62.9883924164, 15.0925357011, 15.0925357008},
 {12.761626995, 12.7616269916, 62.9847225783, 62.9847225786, 12.7616269932, 12.761626993},
 {10.626738228, 10.6267382246, 62.9793388097, 62.97933881, 10.6267382263, 10.6267382261},
 {8.69317519432, 8.69317519065, 62.9720840594, 62.9720840597, 8.69317519232, 8.6931751923},
 {6.96421788614, 6.96421788268, 62.96281659, 62.9628165876, 6.96421788429, 6.96421788429},
 {5.42197032191, 5.42197031861, 62.9515604455, 62.9515604664, 5.42197032016, 5.42197032015},
 {4.18867956307, 4.18867956012, 62.9372793052, 62.937279333, 4.18867956158, 4.18867956157}};

    TestHelper::dumpGolden("correctWaterContentSolution", waterContentSolution);
    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{{0.117508643867, 0.117508643867, 0.080868293801, 0.080868293801, 0.117508643867, 0.117508643867},
 {0.19729680497, 0.19729680497, 0.161029873005, 0.161029873005, 0.19729680497, 0.19729680497},
 {0.26152632582, 0.26152632582, 0.230238760705, 0.230238760705, 0.26152632582, 0.26152632582},
 {0.315643786123, 0.315643786123, 0.289060228318, 0.289060228318, 0.315643786123, 0.315643786123},
 {0.362058093515, 0.362058093515, 0.339333931089, 0.339333931089, 0.362058093515, 0.362058093515},
 {0.402391642418, 0.402391642418, 0.382752993589, 0.382752993589, 0.402391642418, 0.402391642418},
 {0.437895332001, 0.437895332001, 0.420708080501, 0.420708080501, 0.437895332002, 0.437895332002},
 {0.469561960778, 0.469561960778, 0.454317312556, 0.454317312556, 0.469561960779, 0.469561960779},
 {0.498183150342, 0.498183150342, 0.484474973498, 0.484474973498, 0.498183150343, 0.498183150343},
 {0.524389991705, 0.524389991706, 0.511894787524, 0.511894787524, 0.524389991705, 0.524389991706},
 {0.548685052099, 0.548685052099, 0.537145522121, 0.537145522121, 0.548685052098, 0.548685052099},
 {0.571468058787, 0.571468058786, 0.560679824862, 0.560679824862, 0.571468058786, 0.571468058786},
 {0.593056553201, 0.593056553201, 0.582857473568, 0.582857473568, 0.593056553201, 0.593056553201},
 {0.613702484048, 0.613702484048, 0.603964078017, 0.603964078017, 0.613702484048, 0.613702484048},
 {0.638659087722, 0.638659087722, 0.627710208699, 0.627710208699, 0.638659087722, 0.638659087722},
 {0.667336970313, 0.667336970316, 0.654844955951, 0.654844955951, 0.667336970316, 0.667336970316},
 {0.698620077572, 0.698620077575, 0.684827686843, 0.684827686842, 0.698620077575, 0.698620077575},
 {0.80579659564, 0.805796595695, 0.768103940487, 0.768103940482, 0.80579659566, 0.805796595665},
 {0.977614913455, 0.977614913523, 0.912540173995, 0.912540173988, 0.97761491349, 0.977614913496},
 {1.22614269392, 1.22614269403, 1.12885922628, 1.12885922627, 1.22614269398, 1.22614269398},
 {1.54015344733, 1.54015344745, 1.41259125929, 1.41259125928, 1.5401534474, 1.5401534474},
 {1.91142939982, 1.91142939996, 1.75664449417, 1.75664449416, 1.91142939989, 1.9114293999},
 {2.32083188745, 2.32083188759, 2.14557803308, 2.14557803306, 2.32083188753, 2.32083188753},
 {2.83701234569, 2.83701234592, 2.62211679702, 2.62211679699, 2.83701234581, 2.83701234581}};

    TestHelper::dumpGolden("correctTemperatureSolution", temperatureSolution);
    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-6);
}
