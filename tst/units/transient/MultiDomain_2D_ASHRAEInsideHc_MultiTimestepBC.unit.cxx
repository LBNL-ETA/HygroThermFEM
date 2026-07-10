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

    std::vector<std::vector<double>> correctWaterContentSolution{{3.06727324, 3.06727324, 0.00297621939, 0.00297621939, 5.77572618e-06, 5.77572618e-06},
 {4.93776026, 4.93776026, 0.0091609742, 0.0091609742, 2.61145696e-05, 2.61145696e-05},
 {5.96763479, 5.96763479, 0.0170386985, 0.0170386985, 6.86438736e-05, 6.86438736e-05},
 {6.32546516, 6.32546516, 0.0257277192, 0.0257277192, 0.0001389215, 0.0001389215},
 {6.11325712, 6.11325712, 0.0346450564, 0.0346450564, 0.000240183969, 0.000240183969},
 {5.79600933, 5.79600933, 0.0435697906, 0.0435697906, 0.000373984899, 0.000373984899},
 {5.41790473, 5.41790473, 0.0523354102, 0.0523354102, 0.000540341529, 0.000540341529},
 {5.03712762, 5.03712762, 0.0606304276, 0.0606304276, 0.000737372332, 0.000737372332},
 {4.68408218, 4.68408218, 0.0683383251, 0.0683383251, 0.000962208171, 0.000962208171},
 {4.37607283, 4.37607283, 0.0754902469, 0.0754902469, 0.00121170166, 0.00121170166}};

    TestHelper::dumpGolden("correctWaterContentSolution", waterContentSolution);
    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{{3.52432601, 3.52432601, 1.83071368, 1.83071368, 1.37015954, 1.37015954},
 {5.04274486, 5.04274486, 3.43815083, 3.43815083, 2.91790127, 2.91790127},
 {5.98561562, 5.98561562, 4.69222651, 4.69222651, 4.24585081, 4.24585081},
 {6.64621836, 6.64621836, 5.64761452, 5.64761452, 5.2949609, 5.2949609},
 {7.10451025, 7.10451025, 6.3571129, 6.3571129, 6.0898928, 6.0898928},
 {7.39914863, 7.39914863, 6.86236385, 6.86236385, 6.66801674, 6.66801674},
 {7.55802093, 7.55802093, 7.19729511, 7.19729511, 7.06412645, 7.06412645},
 {7.61703879, 7.61703879, 7.3969916, 7.3969916, 7.31323352, 7.31323352},
 {7.59479693, 7.59479693, 7.487948, 7.487948, 7.44397516, 7.44397516},
 {7.50904552, 7.50904552, 7.4923861, 7.4923861, 7.48018701, 7.48018701}};

    TestHelper::dumpGolden("correctTemperatureSolution", temperatureSolution);
    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-6);
}
