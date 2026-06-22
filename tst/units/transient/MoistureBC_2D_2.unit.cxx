#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "TestHelpers.hxx"
#include "TestMaterials.hxx"

TEST(MoistureBC_2D_2, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Simple two elements example with moisture transfer.");

    HygroThermFEM::MultiDomain multiDomain({.performThermal = false, .performMoisture = true});

    std::vector<double> gridXCoordinates{0.15, 0.05, 0.00};

    constexpr HygroThermFEM::State state({
        .temperature = 20.0,
        .humidity = 0.0,
        .pressure = 101325.0,
        .liquidPercent = 1.0
    });

    for(auto val : gridXCoordinates)
    {
        multiDomain.nodes().createNode({.x = val, .y = 0.05, .state = state});
        multiDomain.nodes().createNode({.x = val, .y = 0.00, .state = state});
    }

    const auto & material = multiDomain.materials().createSolidMaterial(TestHelper::CottaerSandstone());

    /// Create elements
    for(size_t i = 1; i <= (multiDomain.nodes().maxIndex() - 2) / 2; ++i)
    {
        const auto node1 = 2u * i + 1u;
        const auto node2 = 2u * i + 2u;
        const auto node3 = 2u * i;
        const auto node4 = 2u * i - 1u;
        multiDomain.createElement({.node1 = node2, .node2 = node3, .node3 = node4, .node4 = node1, .material = material.name()});
    }

    // Create Boundary Conditions
    constexpr auto ambientTemperature = 20.0;
    constexpr auto ambientHumidity = 0.2;
    const auto surfaceTilt{90.0};

    const HygroThermFEM::TARPCoefficients bcCoeff{ambientTemperature, ambientHumidity};

    multiDomain.moisture().createBC_TARPHc(5, 6, bcCoeff, surfaceTilt);

    constexpr auto dTime = 3600;
    constexpr auto nSteps = 4;

    auto humidities = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);
    std::vector<std::vector<double>> solution;

    for(size_t i = 0u; i < nSteps; ++i)
    {
        humidities = multiDomain.moisture().transient(humidities, dTime).value().solution;
        auto waterContent = multiDomain.nodes().properties(HygroThermFEM::Variable::water);
        solution.push_back(waterContent);
    }

    std::vector<std::vector<double>> correctSolution{
      {3.748586e-06, 3.748586e-06, 0.001924969331, 0.001924969331, 0.7428227662, 0.7428227662},
      {1.364027149e-05, 1.364027149e-05, 0.005083314705, 0.005083314705, 1.221651619, 1.221651619},
      {3.11968812e-05, 3.11968812e-05, 0.009029288908, 0.009029288908, 1.530317721, 1.530317721},
      {5.737783803e-05, 5.737783803e-05, 0.01347560965, 0.01347560965, 1.729301858, 1.729301858}};

    TestHelper::expectNear(correctSolution, solution, 1e-9);
}
