#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "SlabCreator.hxx"
#include "TestHelpers.hxx"
#include "TestMaterials.hxx"

using HygroThermFEM::Nodes;

//////////////////////////////////////////////////////////////////////////////////////////////////
// Test case with multi domain where only humidity calculations are performed while temperature is
// kept identical
//////////////////////////////////////////////////////////////////////////////////////////////////

TEST(MultiDomainHumidityOnly_1, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Simple two elements example with moisture and heat transfer.");

    HygroThermFEM::MultiDomain multiDomain({.performThermal = false, .performMoisture = true});

    const auto & material = multiDomain.materials().createSolidMaterial(TestHelper::CottaerSandstone());

    constexpr HygroThermFEM::State state({
        .temperature = 0.0,
        .humidity = 0.6,
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

    std::vector<std::vector<double>> correctWaterContentSolution{{7.0518773, 7.0518773, 7.36650982, 7.36650982, 7.36666651, 7.36666651},
 {6.74538467, 6.74538467, 7.36620042, 7.36620042, 7.36666605, 7.36666605},
 {6.44696994, 6.44696994, 7.36574264, 7.36574265, 7.36666513, 7.36666513},
 {6.15642006, 6.15642006, 7.36514057, 7.36514057, 7.36666361, 7.36666361},
 {5.87352758, 5.87352757, 7.36439816, 7.36439816, 7.36666135, 7.36666135},
 {5.59809053, 5.59809053, 7.36351926, 7.36351926, 7.36665822, 7.36665822},
 {5.3299123, 5.32991229, 7.36250763, 7.36250763, 7.36665408, 7.36665408},
 {5.07444872, 5.07444872, 7.36126302, 7.36126302, 7.36664871, 7.36664871},
 {4.83178951, 4.83178951, 7.35978393, 7.35978393, 7.36664187, 7.36664187},
 {4.60129408, 4.60129407, 7.35808243, 7.35808243, 7.36663334, 7.36663334}};

    TestHelper::dumpGolden("correctWaterContentSolution", waterContentSolution);
    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

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

    TestHelper::dumpGolden("correctTemperatureSolution", temperatureSolution);
    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-6);
}
