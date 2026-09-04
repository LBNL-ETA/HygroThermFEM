#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "TestMaterials.hxx"
#include "TestHelpers.hxx"

TEST(MoistureBC_2D_4, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Simple two elements example with moisture transfer.");

    HygroThermFEM::MultiDomain multiDomain({.performThermal = false, .performMoisture = true});

    std::vector<double> gridXCoordinates{0.15, 0.05, 0.00};

    constexpr HygroThermFEM::State state({
        .temperature = 0.0,
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

    constexpr auto dTime = 36000;
    constexpr auto nSteps = 4;

    auto humidities = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);
    std::vector<std::vector<double>> solution;

    for(size_t i = 0u; i < nSteps; ++i)
    {
        humidities = multiDomain.moisture().transient(humidities, dTime).value().solution;
        auto waterContent = multiDomain.nodes().properties(HygroThermFEM::Variable::water);
        solution.push_back(waterContent);
    }

    std::vector<std::vector<double>> correctSolution{{0.000150023209022, 0.000150023209022, 0.0310126698632, 0.0310126698632, 4.83136883365, 4.83136883365},
 {0.000498024227103, 0.000498024227104, 0.0720884946388, 0.0720884946388, 7.5316035421, 7.5316035421},
 {0.00106025188291, 0.00106025188291, 0.116721239685, 0.116721239685, 8.78226358422, 8.78226358422},
 {0.00184128368075, 0.00184128368075, 0.162514467027, 0.162514467027, 9.34009200043, 9.34009200043}};

    TestHelper::dumpGolden("correctSolution", solution);
    EXPECT_EQ(solution.size(), correctSolution.size());

    std::cout.precision(10);
    for(auto i = 0u; i < correctSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctSolution[i][j], solution[i][j], 1e-8);
        }
    }
}
