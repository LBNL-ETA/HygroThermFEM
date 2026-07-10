#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "SlabCreator.hxx"
#include "TestHelpers.hxx"
#include "TestMaterials.hxx"

TEST(MultiDomain_2D_ASHRAEOutsideHc_MultiTimestepBC, TestExample_1)
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

    // Variable boundary conditions (temperature, humidity and wind speed) over ten timesteps.
    const std::vector<HygroThermFEM::ASHRAEOutsideCoefficients> bcCoeff{{20.0, 0.6, 3},
                                                                       {20.0, 0.5, 3},
                                                                       {20.0, 0.4, 3},
                                                                       {20.0, 0.3, 4},
                                                                       {20.0, 0.2, 4.2},
                                                                       {18.0, 0.2, 4.6},
                                                                       {16.0, 0.2, 5},
                                                                       {14.0, 0.2, 5.3},
                                                                       {12.0, 0.2, 5.5},
                                                                       {10.0, 0.2, 5.9}};

    multiDomain.createBC_ASHRAEOutsideHc(1, 2, bcCoeff);

    constexpr auto dTime = 3600;
    constexpr auto nSteps = 10;

    const auto results = multiDomain.transientMultiStep(dTime, nSteps);

    std::vector<std::vector<double>> correctWaterContentSolution{{11.6503249288, 11.6503249288, 0.00809627974437, 0.00809627974437, 1.57118440916e-05, 1.57118440916e-05},
 {12.2914824044, 12.2914824044, 0.0266132866794, 0.0266132866794, 9.48895996264e-05, 9.48895996264e-05},
 {10.3481974351, 10.3481974351, 0.0445953528111, 0.0445953528111, 0.00026097065829, 0.00026097065829},
 {7.09408719031, 7.09408719031, 0.059513970143, 0.059513970143, 0.000512818952412, 0.000512818952412},
 {4.16832398649, 4.16832398649, 0.0700345530101, 0.0700345530101, 0.00083897673322, 0.00083897673322},
 {2.9615545341, 2.9615545341, 0.0778621006338, 0.0778621006338, 0.00123128054013, 0.00123128054013},
 {2.37749769593, 2.37749769593, 0.0842430230251, 0.0842430230251, 0.00167919672954, 0.00167919672954},
 {2.06531468332, 2.06531468332, 0.0896128585301, 0.0896128585301, 0.00216387970349, 0.00216387970349},
 {1.88087104505, 1.88087104505, 0.0941894850772, 0.0941894850772, 0.00266394563126, 0.00266394563126},
 {1.76119878493, 1.76119878493, 0.0981058766481, 0.0981058766481, 0.00316013802822, 0.00316013802822}};

    TestHelper::dumpGolden("correctWaterContentSolution", results.moisture.values);
    TestHelper::expectNear(correctWaterContentSolution, results.moisture.values, 1e-6);

    const std::vector<double> correctHumidityError{3.53465703447e-07, 1.53019311717e-08, 1.3152823237e-07, 3.12934512626e-07, 6.02696193286e-07, 5.00279517585e-07, 2.0238109354e-07, 2.39875799451e-09, 1.26541530036e-07, 2.12035373267e-07};
    TestHelper::dumpGolden("correctHumidityError", results.moisture.errors);
    TestHelper::expectNear(correctHumidityError, results.moisture.errors, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{{11.3429351181, 11.3429351181, 5.88913370454, 5.88913370454, 4.40760150594, 4.40760150594},
 {12.393187695, 12.393187695, 9.07902701395, 9.07902701395, 7.9038432106, 7.9038432106},
 {12.95244842, 12.95244842, 10.9445028075, 10.9445028075, 10.1795997262, 10.1795997262},
 {13.8958700592, 13.8958700592, 12.3842802356, 12.3842802356, 11.8297155018, 11.8297155018},
 {14.8514011014, 14.8514011014, 13.5984964197, 13.5984964197, 13.1536183522, 13.1536183522},
 {15.3558411809, 15.3558411809, 14.4572181064, 14.4572181064, 14.1293928727, 14.1293928727},
 {15.1406164238, 15.1406164238, 14.772655107, 14.772655107, 14.6109792212, 14.6109792212},
 {14.3889967174, 14.3889967174, 14.5549398306, 14.5549398306, 14.5692114113, 14.5692114113},
 {13.2672826231, 13.2672826231, 13.8901479858, 13.8901479858, 14.0611706841, 14.0611706841},
 {11.8559048518, 11.8559048518, 12.8576108348, 12.8576108348, 13.1605906944, 13.1605906944}};

    TestHelper::dumpGolden("correctTemperatureSolution", results.temperature.values);
    TestHelper::expectNear(correctTemperatureSolution, results.temperature.values, 1e-6);

    const std::vector<double> correctTemperatureError{3.59712485754e-07, 1.77468152291e-07, 1.07427668275e-07, 6.91044690603e-08, 5.97949215248e-08, 4.35083369649e-08, 1.60375833568e-08, 1.29882302118e-08, 4.12552915975e-08, 7.01695531172e-08};
    TestHelper::dumpGolden("correctTemperatureError", results.temperature.errors);
    TestHelper::expectNear(correctTemperatureError, results.temperature.errors, 1e-6);
}
