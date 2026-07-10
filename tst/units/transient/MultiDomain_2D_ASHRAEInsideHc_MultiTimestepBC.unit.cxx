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
 {4.93790395548, 4.93790395548, 0.00916035982241, 0.00916035982241, 2.61144800634e-05, 2.61144800634e-05},
 {5.96811172046, 5.96811172046, 0.017036247443, 0.017036247443, 6.8643967885e-05, 6.8643967885e-05},
 {6.32631204075, 6.32631204074, 0.0257228705522, 0.0257228705522, 0.000138923643458, 0.000138923643458},
 {6.11441398264, 6.11441398263, 0.0346378104782, 0.0346378104782, 0.000240192390232, 0.000240192390232},
 {5.79740304434, 5.79740304434, 0.0435604532453, 0.0435604532453, 0.00037400618077, 0.00037400618077},
 {5.41945512766, 5.41945512766, 0.0523243725787, 0.0523243725787, 0.000540384346595, 0.000540384346595},
 {5.03862302431, 5.03862302431, 0.060619250703, 0.060619250703, 0.000737450218007, 0.000737450218007},
 {4.68547817621, 4.68547817621, 0.0683272718075, 0.0683272718075, 0.00096233304166, 0.00096233304166},
 {4.3773244109, 4.37732441091, 0.0754794908615, 0.0754794908615, 0.00121188361373, 0.00121188361373}};

    TestHelper::dumpGolden("correctWaterContentSolution", waterContentSolution);
    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{{3.52247587826, 3.52247587826, 1.83156898532, 1.83156898532, 1.37080222256, 1.37080222256},
 {5.03871008649, 5.03871008649, 3.44019584309, 3.44019584309, 2.91960244167, 2.91960244167},
 {5.98044757642, 5.98044757642, 4.69525707906, 4.69525707906, 4.24856561424, 4.24856561424},
 {6.64067720657, 6.64067720657, 5.65141883882, 5.65141883882, 5.29852195532, 5.29852195532},
 {7.09907304731, 7.09907304731, 6.36149500508, 6.36149500508, 6.09411273148, 6.09411273148},
 {7.39399337133, 7.39399337133, 6.86719170484, 6.86719170484, 6.67275026548, 6.67275026548},
 {7.55323861118, 7.55323861118, 7.20245137442, 7.20245137442, 7.0692491947, 7.0692491947},
 {7.61272385976, 7.61272385976, 7.40228067112, 7.40228067112, 7.3185669833, 7.3185669833},
 {7.59104433956, 7.59104433956, 7.49325586326, 7.49325586326, 7.44938786004, 7.44938786004},
 {7.50580442941, 7.50580442941, 7.4976473701, 7.4976473701, 7.48559554907, 7.48559554907}};

    TestHelper::dumpGolden("correctTemperatureSolution", temperatureSolution);
    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-6);
}
