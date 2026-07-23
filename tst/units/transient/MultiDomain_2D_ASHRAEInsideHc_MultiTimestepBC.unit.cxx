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

    std::vector<std::vector<double>> correctWaterContentSolution{{3.06727324303, 3.06727324303, 0.00297621939188, 0.00297621939188, 5.77572618236e-06, 5.77572618236e-06},
 {4.93793660598, 4.93793660598, 0.00915791337621, 0.00915791337621, 2.6108428395e-05, 2.6108428395e-05},
 {5.96824882213, 5.96824882213, 0.0170306469504, 0.0170306469504, 6.86215430778e-05, 6.86215430778e-05},
 {6.32659128549, 6.32659128548, 0.0257145127471, 0.0257145127471, 0.000138873907493, 0.000138873907493},
 {6.11483428104, 6.11483428103, 0.0346272453414, 0.0346272453414, 0.000240105206499, 0.000240105206499},
 {5.79794954556, 5.79794954556, 0.0435481583165, 0.0435481583166, 0.000373872811844, 0.000373872811844},
 {5.42010014653, 5.42010014653, 0.0523107337855, 0.0523107337855, 0.000540197811912, 0.000540197811912},
 {5.03927324378, 5.03927324377, 0.0606050924496, 0.0606050924496, 0.000737207122281, 0.000737207122281},
 {4.68610617448, 4.68610617446, 0.0683128003403, 0.0683128003403, 0.000962031259778, 0.000962031259778},
 {4.37789637225, 4.37789637223, 0.0754648167731, 0.0754648167732, 0.00121152202075, 0.00121152202075}};

    TestHelper::dumpGolden("correctWaterContentSolution", waterContentSolution);
    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{{3.52210584038, 3.52210584038, 1.83155828414, 1.83155828414, 1.37079434535, 1.37079434535},
 {5.03751417924, 5.03751417924, 3.43992820928, 3.43992820928, 2.91940089897, 2.91940089897},
 {5.97856078829, 5.97856078829, 4.69474064842, 4.69474064842, 4.24813054616, 4.24813054616},
 {6.63833742218, 6.63833742218, 5.65073873416, 5.65073873416, 5.29790785886, 5.29790785886},
 {7.09649998859, 7.09649998859, 6.36073549146, 6.36073549146, 6.09339698977, 6.09339698977},
 {7.39133213259, 7.39133213259, 6.86642025418, 6.86642025418, 6.67200318989, 6.67200318989},
 {7.55059555802, 7.55059555802, 7.20171360444, 7.20171360444, 7.06852269367, 7.06852269367},
 {7.61016502538, 7.61016502538, 7.4015787102, 7.4015787102, 7.31787551526, 7.31787551526},
 {7.58863667603, 7.58863667603, 7.49259573415, 7.49259573415, 7.4487391956, 7.4487391956},
 {7.50355423501, 7.50355423501, 7.49702723296, 7.49702723296, 7.48498981447, 7.48498981447}};

    TestHelper::dumpGolden("correctTemperatureSolution", temperatureSolution);
    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-6);
}
