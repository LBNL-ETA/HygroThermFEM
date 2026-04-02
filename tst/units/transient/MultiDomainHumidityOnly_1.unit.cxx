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

    std::vector<std::vector<double>> correctWaterContentSolution{
      {7.055972, 7.055972, 7.366497, 7.366497, 7.366666, 7.366666},
      {6.748959, 6.748959, 7.366159, 7.366159, 7.366666, 7.366666},
      {6.444629, 6.444629, 7.365653, 7.365653, 7.366665, 7.366665},
      {6.143056, 6.143056, 7.364982, 7.364982, 7.366663, 7.366663},
      {5.844321, 5.844321, 7.364145, 7.364145, 7.366660, 7.366660},
      {5.548503, 5.548503, 7.363145, 7.363145, 7.366656, 7.366656},
      {5.277271, 5.277271, 7.361982, 7.361982, 7.366651, 7.366651},
      {5.128666, 5.128666, 7.360658, 7.360658, 7.366645, 7.366645},
      {4.981688, 4.981688, 7.359174, 7.359174, 7.366637, 7.366637},
      {4.836384, 4.836384, 7.357531, 7.357531, 7.366627, 7.366627}};

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

    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-6);
}
