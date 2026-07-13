#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "SlabCreator.hxx"
#include "TestHelpers.hxx"
#include "TestMaterials.hxx"

TEST(MultiDomain_2D_VariableTARPHc_MultiTimestepBC, TestExample_1)
{
    HygroThermFEM::MultiDomain multiDomain;

    const auto & material =
      multiDomain.materials().createSolidMaterial(TestHelper::CottaerSandstone());

    constexpr HygroThermFEM::State state(
      {.temperature = 0.0, .humidity = 0.0, .pressure = 101325.0, .liquidPercent = 1.0});

    TestHelper::SlabBuilder(multiDomain)
      .gridXCoordinates({0, 0.05, 0.1})
      .height(0.05)
      .material(material.name())
      .state(state)
      .startCorner(TestHelper::StartCorner::BottomRight)
      .direction(TestHelper::Direction::CounterClockwise)
      .build();

    /// Create Boundary Conditions

    // Variable boundary conditions (temperature and humidity) over ten timesteps.
    const std::vector<HygroThermFEM::TARPCoefficients> bcCoeff{{20.0, 0.6},
                                                               {20.0, 0.5},
                                                               {20.0, 0.4},
                                                               {20.0, 0.3},
                                                               {20.0, 0.2},
                                                               {18.0, 0.2},
                                                               {16.0, 0.2},
                                                               {14.0, 0.2},
                                                               {12.0, 0.2},
                                                               {10.0, 0.2}};

    const auto surfaceTilt{90.0};

    multiDomain.createBC_TARPHc(1, 2, bcCoeff, surfaceTilt);

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

    std::vector<std::vector<double>> correctWaterContentSolution{{3.46405340748, 3.46405340748, 0.00336122089849, 0.00336122089849, 6.5228697861e-06, 6.5228697861e-06},
 {5.44168826425, 5.44168826425, 0.0103050942738, 0.0103050942738, 2.97894304346e-05, 2.97894304346e-05},
 {6.4844989023, 6.48449890232, 0.0188458558194, 0.0188458558194, 7.80822450734e-05, 7.80822450734e-05},
 {6.77751527699, 6.777515277, 0.0282097989918, 0.0282097989918, 0.000157597538806, 0.000157597538806},
 {6.45419319757, 6.45419319759, 0.0377258887932, 0.0377258887932, 0.000271629841029, 0.000271629841029},
 {6.0234160626, 6.02341606262, 0.047146641217, 0.047146641217, 0.000421424698162, 0.000421424698162},
 {5.51676555315, 5.51676555316, 0.0562920213711, 0.056292021371, 0.000606597120289, 0.000606597120289},
 {4.99917084596, 4.99917084597, 0.06480262794, 0.0648026279399, 0.000824637822601, 0.000824637822601},
 {4.51168825814, 4.51168825815, 0.0724731188741, 0.0724731188741, 0.00107161281774, 0.00107161281774},
 {4.06151143843, 4.06151143845, 0.0793182468618, 0.0793182468617, 0.00134312554706, 0.00134312554706}};

    TestHelper::dumpGolden("correctWaterContentSolution", waterContentSolution);
    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{{3.96182969328, 3.96182969328, 2.06012547793, 2.06012547793, 1.54186093011, 1.54186093011},
 {5.53518426592, 5.53518426592, 3.801065778, 3.801065778, 3.23272280168, 3.23272280168},
 {6.49388112119, 6.49388112119, 5.13001850517, 5.13001850517, 4.65272956273, 4.65272956273},
 {7.14577030769, 7.14577030769, 6.11986207202, 6.11986207202, 5.75080216732, 5.75080216732},
 {7.58694462483, 7.58694462483, 6.83919485241, 6.83919485241, 6.56543027078, 6.56543027078},
 {7.88322978784, 7.88322978784, 7.35137009968, 7.35137009968, 7.15371076552, 7.15371076552},
 {8.05437923841, 8.05437923841, 7.69631801459, 7.69631801459, 7.55989140028, 7.55989140028},
 {8.12131013927, 8.12131013927, 7.90455082595, 7.90455082595, 7.81793863634, 7.81793863634},
 {8.09609643627, 8.09609643627, 7.99755862657, 7.99755862657, 7.95248040044, 7.95248040044},
 {7.97739614247, 7.97739614247, 7.98552937418, 7.98552937418, 7.9773370164, 7.9773370164}};

    TestHelper::dumpGolden("correctTemperatureSolution", temperatureSolution);
    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-6);
}
