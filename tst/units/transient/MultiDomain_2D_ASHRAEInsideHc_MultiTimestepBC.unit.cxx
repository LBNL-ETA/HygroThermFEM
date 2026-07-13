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
 {4.93793269928, 4.93793269928, 0.00916023422243, 0.00916023422243, 2.61141850765e-05, 2.61141850765e-05},
 {5.96824056116, 5.96824056116, 0.01703555286, 0.01703555286, 6.86411141907e-05, 6.86411141907e-05},
 {6.32658048884, 6.32658048884, 0.0257212387671, 0.0257212387671, 0.000138913722453, 0.000138913722453},
 {6.11482252847, 6.11482252847, 0.0346350605681, 0.0346350605681, 0.000240169512778, 0.000240169512778},
 {5.79793779165, 5.79793779165, 0.04355656234, 0.04355656234, 0.000373964139667, 0.000373964139667},
 {5.4200889059, 5.4200889059, 0.0523194142297, 0.0523194142297, 0.000540317481546, 0.000540317481546},
 {5.03926359442, 5.03926359442, 0.0606138604424, 0.0606138604424, 0.000737355635346, 0.000737355635346},
 {4.68609788633, 4.68609788632, 0.0683215705341, 0.0683215705341, 0.00096220869543, 0.00096220869543},
 {4.37788918802, 4.37788918802, 0.0754735558992, 0.0754735558992, 0.00121172823326, 0.00121172823326}};

    TestHelper::dumpGolden("correctWaterContentSolution", waterContentSolution);
    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{{3.52210584038, 3.52210584038, 1.83155828414, 1.83155828414, 1.37079434535, 1.37079434535},
 {5.03751457995, 5.03751457995, 3.43992840762, 3.43992840762, 2.91940104863, 2.91940104863},
 {5.97856147923, 5.97856147923, 4.69474107974, 4.69474107974, 4.24813091015, 4.24813091015},
 {6.63833849746, 6.63833849746, 5.65073947038, 5.65073947038, 5.29790850754, 5.29790850754},
 {7.09650142826, 7.09650142826, 6.36073655958, 6.36073655958, 6.09339796072, 6.09339796072},
 {7.39133386419, 7.39133386419, 6.86642163266, 6.86642163266, 6.67200447593, 6.67200447593},
 {7.55059751487, 7.55059751487, 7.20171525108, 7.20171525108, 7.06852426088, 7.06852426088},
 {7.61016760705, 7.61016760705, 7.4015808092, 7.4015808092, 7.31787749253, 7.31787749253},
 {7.58863953165, 7.58863953165, 7.49259819011, 7.49259819011, 7.44874154362, 7.44874154362},
 {7.5035572458, 7.5035572458, 7.49702994523, 7.49702994523, 7.48499244782, 7.48499244782}};

    TestHelper::dumpGolden("correctTemperatureSolution", temperatureSolution);
    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-6);
}
