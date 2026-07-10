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

    std::vector<std::vector<double>> correctWaterContentSolution{{60.3402862, 60.3402862, 62.9999423, 62.9999423, 60.3402862, 60.3402862},
 {57.666585, 57.666585, 62.9998775, 62.9998775, 57.666585, 57.666585},
 {54.9859115, 54.9859115, 62.9997528, 62.9997528, 54.9859115, 54.9859115},
 {52.3011719, 52.3011719, 62.9995609, 62.9995609, 52.3011719, 52.3011719},
 {49.6142577, 49.6142577, 62.9993018, 62.9993018, 49.6142577, 49.6142577},
 {46.9266089, 46.9266089, 62.9989763, 62.9989763, 46.9266089, 46.9266089},
 {44.2393603, 44.2393603, 62.9985852, 62.9985852, 44.2393603, 44.2393603},
 {41.5534106, 41.5534106, 62.9981292, 62.9981292, 41.5534106, 41.5534106},
 {38.8694775, 38.8694775, 62.997609, 62.997609, 38.8694775, 38.8694775},
 {36.1881303, 36.1881303, 62.9970249, 62.9970249, 36.1881303, 36.1881303},
 {33.5098188, 33.5098188, 62.9963774, 62.9963774, 33.5098188, 33.5098188},
 {30.834898, 30.834898, 62.9956668, 62.9956668, 30.834898, 30.834898},
 {28.1636473, 28.1636473, 62.9948932, 62.9948932, 28.1636473, 28.1636473},
 {25.4962863, 25.4962863, 62.9940569, 62.9940569, 25.4962863, 25.4962863},
 {22.8417954, 22.8417954, 62.9930908, 62.9930908, 22.8417954, 22.8417954},
 {20.2013316, 20.2013316, 62.9919817, 62.9919817, 20.2013316, 20.2013316},
 {17.5741644, 17.5741644, 62.9907298, 62.9907298, 17.5741644, 17.5741644},
 {15.0869902, 15.0869902, 62.9883646, 62.9883646, 15.0869902, 15.0869902},
 {12.7549432, 12.7549432, 62.9846861, 62.9846861, 12.7549432, 12.7549432},
 {10.6185279, 10.6185279, 62.9792842, 62.9792842, 10.6185279, 10.6185279},
 {8.68227417, 8.68227417, 62.9719965, 62.9719965, 8.68227417, 8.68227417},
 {6.94929463, 6.94929463, 62.9626661, 62.9626661, 6.94929463, 6.94929463},
 {5.40103529, 5.40103529, 62.9512974, 62.9512974, 5.40103529, 5.40103529},
 {4.16380538, 4.16380538, 62.9367898, 62.9367898, 4.16380538, 4.16380538}};

    TestHelper::dumpGolden("correctWaterContentSolution", waterContentSolution);
    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{{0.117508702, 0.117508702, 0.0808479344, 0.0808479344, 0.117508702, 0.117508702},
 {0.197299365, 0.197299365, 0.160994092, 0.160994092, 0.197299365, 0.197299365},
 {0.261548937, 0.261548937, 0.230223542, 0.230223542, 0.261548937, 0.261548937},
 {0.31572875, 0.31572875, 0.289116068, 0.289116068, 0.31572875, 0.31572875},
 {0.362264045, 0.362264045, 0.339525321, 0.339525321, 0.362264045, 0.362264045},
 {0.402788526, 0.402788526, 0.383156289, 0.383156289, 0.402788526, 0.402788526},
 {0.43856047, 0.43856047, 0.421408008, 0.421408008, 0.43856047, 0.43856047},
 {0.470537486, 0.470537486, 0.455362645, 0.455362645, 0.470537486, 0.470537486},
 {0.499507909, 0.499507909, 0.485909166, 0.485909166, 0.499507909, 0.499507909},
 {0.526102525, 0.526102525, 0.513760481, 0.513760481, 0.526102525, 0.526102525},
 {0.550823277, 0.550823277, 0.539484647, 0.539484647, 0.550823277, 0.550823277},
 {0.574069299, 0.574069299, 0.563533686, 0.563533686, 0.574069299, 0.574069299},
 {0.596157676, 0.596157676, 0.586266892, 0.586266892, 0.596157676, 0.596157676},
 {0.617440089, 0.617440089, 0.608073041, 0.608073041, 0.617440089, 0.617440089},
 {0.643219885, 0.643219885, 0.632748324, 0.632748324, 0.643219885, 0.643219885},
 {0.672877921, 0.672877921, 0.661000737, 0.661000737, 0.672877921, 0.672877921},
 {0.705292978, 0.705292978, 0.692272472, 0.692272472, 0.705292978, 0.705292978},
 {0.816059642, 0.816059642, 0.779726311, 0.779726311, 0.816059642, 0.816059642},
 {0.994909548, 0.994909548, 0.932218952, 0.932218952, 0.994909548, 0.994909548},
 {1.25673593, 1.25673593, 1.16335058, 1.16335058, 1.25673593, 1.25673593},
 {1.59293281, 1.59293281, 1.47149096, 1.47149096, 1.59293281, 1.59293281},
 {1.99895303, 1.99895303, 1.85332696, 1.85332696, 1.99895303, 1.99895303},
 {2.45824914, 2.45824914, 2.29626844, 2.29626844, 2.45824914, 2.45824914},
 {3.05379745, 3.05379745, 2.8572566, 2.8572566, 3.05379745, 3.05379745}};

    TestHelper::dumpGolden("correctTemperatureSolution", temperatureSolution);
    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-6);
}
