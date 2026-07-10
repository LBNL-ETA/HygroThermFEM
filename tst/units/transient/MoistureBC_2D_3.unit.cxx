#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "SlabCreator.hxx"
#include "TestHelpers.hxx"
#include "TestMaterials.hxx"

TEST(MoistureBC_2D_3, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Simple two elements example with moisture transfer.");

    HygroThermFEM::MultiDomain multiDomain({.performThermal = false, .performMoisture = true});

    constexpr HygroThermFEM::State state({
        .temperature = 20.0,
        .humidity = 0.0,
        .pressure = 0.0,
        .liquidPercent = 0
    });

    const auto & material = multiDomain.materials().createSolidMaterial(TestHelper::CottaerSandstone());

    TestHelper::SlabBuilder(multiDomain)
        .gridXCoordinates({0, 0.05, 0.1, 0.15})
        .height(0.05)
        .material(material.name())
        .state(state)
        .startCorner(TestHelper::StartCorner::BottomRight)
        .direction(TestHelper::Direction::CounterClockwise)
        .build();

    // Create Boundary Conditions
    constexpr auto ambientTemperature = 20.0;
    constexpr auto ambientHumidity = 0.5;
    const auto surfaceTilt{90.0};

    const HygroThermFEM::TARPCoefficients bcCoeff{ambientTemperature, ambientHumidity};

    multiDomain.moisture().createBC_TARPHc(1, 2, bcCoeff, surfaceTilt);

    constexpr auto dTime = 36000;
    constexpr auto nSteps = 4;

    auto humidities = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);
    std::vector<std::vector<double>> solution;

    for(unsigned i = 0; i < nSteps; ++i)
    {
        humidities = multiDomain.moisture().transient(humidities, dTime).value().solution;
        auto waterContent = multiDomain.nodes().properties(HygroThermFEM::Variable::water);
        solution.push_back(waterContent);
    }

    std::vector<std::vector<double>> correctSolution{{4.42476715, 4.42476715, 0.160377873, 0.160377873, 0.00582059189, 0.00582059189, 0.000421386146, 0.000421386146},
 {5.10613707, 5.10613707, 0.334234718, 0.334234718, 0.0175579543, 0.0175579543, 0.00166200094, 0.00166200094},
 {5.2127262, 5.2127262, 0.499978215, 0.499978215, 0.0345312811, 0.0345312811, 0.00404159731, 0.00404159731},
 {5.23094899, 5.23094899, 0.655160458, 0.655160458, 0.0560293173, 0.0560293173, 0.00780528758, 0.00780528758}};

    TestHelper::dumpGolden("correctSolution", solution);
    TestHelper::expectNear(correctSolution, solution, 1e-6);
}
