#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "SlabCreator.hxx"
#include "TestHelpers.hxx"
#include "TestMaterials.hxx"

TEST(MultiDomain_2D_ASHRAEInsideHc_MultiTimestepBC, TestExample_1)
{
    HygroThermFEM::MultiDomain multiDomain;

    constexpr HygroThermFEM::State state({
        .temperature = 0.0,
        .humidity = 0.0,
        .pressure = 101325.0,
        .liquidPercent = 1.0
    });

    const auto & material = multiDomain.materials().createSolidMaterial(TestHelper::CottaerSandstone());

    TestHelper::SlabBuilder(multiDomain)
        .gridXCoordinates({0, 0.05, 0.1})
        .height(0.05)
        .material(material.name())
        .state(state)
        .startCorner(TestHelper::StartCorner::BottomRight)
        .direction(TestHelper::Direction::CounterClockwise)
        .build();

    /// Create Boundary Conditions

    // Variable boundary conditions (temperature, pressure and humidity) over ten timesteps.
    const std::vector<HygroThermFEM::ASHRAEInsideCoefficients> bcCoeff{{20.0, 0.6, 101325.0},
                                                                       {20.0, 0.5, 101325.0},
                                                                       {20.0, 0.4, 101325.0},
                                                                       {20.0, 0.3, 101325.0},
                                                                       {20.0, 0.2, 101325.0},
                                                                       {18.0, 0.2, 101325.0},
                                                                       {16.0, 0.2, 101325.0},
                                                                       {14.0, 0.2, 101325.0},
                                                                       {12.0, 0.2, 101325.0},
                                                                       {10.0, 0.2, 101325.0}};

    constexpr auto surfaceTilt{90.0};    // degrees
    constexpr auto surfaceHeight{1.0};   // meters

    multiDomain.createBC_ASHRAEInsideHc(1, 2, bcCoeff, surfaceHeight, surfaceTilt);

    constexpr auto dTime = 3600;
    constexpr auto nSteps = 10;

    auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
    auto humidities = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);
    std::vector<std::vector<double>> temperatureSolution;
    std::vector<std::vector<double>> waterContentSolution;
    size_t timestepIndex{0u};

    for(auto i = 0; i < nSteps; ++i)
    {
        auto aSolution = multiDomain.transient(temperatures, humidities, dTime, timestepIndex);
        temperatureSolution.push_back(aSolution.temperature);
        waterContentSolution.push_back(aSolution.waterContent);
        temperatures = aSolution.temperature;
        humidities = aSolution.humidity;
        ++timestepIndex;
    }

    std::vector<std::vector<double>> correctWaterContentSolution{{3.0672732428, 3.0672732428, 0.0029762194552, 0.0029762194552, 5.77572630546e-06, 5.77572630546e-06},
 {4.93793656519, 4.93793656519, 0.00915793758677, 0.00915793758688, 2.61084823195e-05, 2.61084823194e-05},
 {5.96824710476, 5.9682471048, 0.0170306432506, 0.0170306432508, 6.86215876293e-05, 6.8621587629e-05},
 {6.32658961875, 6.32658961883, 0.0257145706875, 0.0257145706877, 0.000138874111029, 0.000138874111029},
 {6.11483286207, 6.11483286218, 0.0346272535197, 0.0346272535199, 0.000240105434418, 0.000240105434417},
 {5.79794836917, 5.79794836938, 0.0435481037767, 0.0435481037767, 0.000373872872303, 0.000373872872303},
 {5.42009917782, 5.42009917799, 0.0523106208245, 0.0523106208246, 0.000540197512897, 0.000540197512896},
 {5.03927150424, 5.03927150445, 0.0606049763327, 0.0606049763326, 0.000737206446169, 0.000737206446169},
 {4.68610474671, 4.6861047468, 0.0683126500674, 0.0683126500675, 0.000962030090885, 0.000962030090884},
 {4.3778951683, 4.37789516838, 0.0754646449656, 0.0754646449662, 0.00121152028704, 0.00121152028704}};

    TestHelper::dumpGolden("correctWaterContentSolution", waterContentSolution);
    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{{3.52210584041, 3.52210584041, 1.83155828415, 1.83155828415, 1.37079434536, 1.37079434536},
 {5.03751418349, 5.03751418349, 3.43992821137, 3.43992821137, 2.91940090055, 2.91940090055},
 {5.97856089022, 5.97856089022, 4.6947407019, 4.6947407019, 4.24813058658, 4.24813058658},
 {6.63833758072, 6.63833758072, 5.65073883984, 5.65073883984, 5.29790794818, 5.29790794818},
 {7.09650018045, 7.09650018045, 6.36073563924, 6.36073563924, 6.09339712285, 6.09339712285},
 {7.39133234288, 7.39133234287, 6.86642043215, 6.86642043215, 6.6720033565, 6.6720033565},
 {7.55059577823, 7.55059577823, 7.20171380267, 7.20171380267, 7.06852288379, 7.06852288379},
 {7.61016539254, 7.61016539253, 7.40157899361, 7.40157899361, 7.31787577505, 7.31787577505},
 {7.58863709663, 7.58863709663, 7.49259608459, 7.4925960846, 7.44873952302, 7.44873952302},
 {7.50355468249, 7.50355468248, 7.49702762989, 7.49702762989, 7.48499019367, 7.48499019367}};

    TestHelper::dumpGolden("correctTemperatureSolution", temperatureSolution);
    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-6);
}
