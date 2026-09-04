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

    std::vector<std::vector<double>> correctWaterContentSolution{{7.05187730268, 7.05187730261, 7.36650981763, 7.36650981763, 7.36666651036, 7.36666651036},
 {6.74538466728, 6.74538466731, 7.36620041949, 7.36620041949, 7.36666604589, 7.36666604589},
 {6.44696992903, 6.44696992907, 7.36574265259, 7.36574265267, 7.36666512569, 7.36666512569},
 {6.15642003082, 6.15642003089, 7.36514058539, 7.36514058549, 7.36666360644, 7.36666360644},
 {5.87352753669, 5.87352753669, 7.3643981779, 7.364398178, 7.36666134887, 7.36666134887},
 {5.598090483, 5.59809048304, 7.36351928504, 7.36351928511, 7.3666582177, 7.3666582177},
 {5.32991223398, 5.32991223403, 7.36250765937, 7.36250765941, 7.36665408153, 7.36665408153},
 {5.07444698052, 5.07444698056, 7.36126304644, 7.36126304648, 7.36664870919, 7.36664870919},
 {4.8317878482, 4.83178784822, 7.3597839647, 7.35978396479, 7.36664186825, 7.36664186825},
 {4.6012924849, 4.601292485, 7.35808246339, 7.35808246354, 7.36663333852, 7.36663333852}};

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
