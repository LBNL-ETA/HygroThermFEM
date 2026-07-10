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

    std::vector<std::vector<double>> correctWaterContentSolution{{11.6503249, 11.6503249, 0.00809627974, 0.00809627974, 1.57118441e-05, 1.57118441e-05},
 {12.2892734, 12.2892734, 0.0266161914, 0.0266161914, 9.48810705e-05, 9.48810705e-05},
 {10.3434035, 10.3434035, 0.0446017319, 0.0446017319, 0.000260886661, 0.000260886661},
 {7.08947684, 7.08947684, 0.0595201927, 0.0595201927, 0.000512560979, 0.000512560979},
 {4.16636225, 4.16636225, 0.0700381894, 0.0700381894, 0.000838460119, 0.000838460119},
 {2.9608141, 2.9608141, 0.0778646869, 0.0778646869, 0.00123047876, 0.00123047876},
 {2.37731035, 2.37731035, 0.0842452322, 0.0842452322, 0.00167810457, 0.00167810457},
 {2.06533541, 2.06533541, 0.0896150143, 0.0896150143, 0.00216250617, 0.00216250617},
 {1.88096271, 1.88096271, 0.0941917077, 0.0941917077, 0.00266231225, 0.00266231225},
 {1.76131099, 1.76131099, 0.0981082029, 0.0981082029, 0.00315827473, 0.00315827473}};

    TestHelper::dumpGolden("correctWaterContentSolution", results.moisture.values);
    TestHelper::expectNear(correctWaterContentSolution, results.moisture.values, 1e-6);

    const std::vector<double> correctHumidityError{3.53193856e-07, 1.55744273e-08, 1.31709146e-07, 3.12942675e-07, 6.02708777e-07, 4.9965971e-07, 2.02011885e-07, 2.56150646e-09, 1.26593354e-07, 2.12032261e-07};
    TestHelper::dumpGolden("correctHumidityError", results.moisture.errors);
    TestHelper::expectNear(correctHumidityError, results.moisture.errors, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{{11.3464588, 11.3464588, 5.88607603, 5.88607603, 4.40530613, 4.40530613},
 {12.4004974, 12.4004974, 9.07039158, 9.07039158, 7.89676806, 7.89676806},
 {12.9583876, 12.9583876, 10.9328425, 10.9328425, 10.1690203, 10.1690203},
 {13.8995979, 13.8995979, 12.3717383, 12.3717383, 11.8175574, 11.8175574},
 {14.8532396, 14.8532396, 13.5870685, 13.5870685, 13.1418646, 13.1418646},
 {15.3559791, 15.3559791, 14.4469537, 14.4469537, 14.118583, 14.118583},
 {15.1399677, 15.1399677, 14.7633711, 14.7633711, 14.6011165, 14.6011165},
 {14.3880179, 14.3880179, 14.5465511, 14.5465511, 14.5602411, 14.5602411},
 {13.2661948, 13.2661948, 13.8826056, 13.8826056, 14.0530515, 14.0530515},
 {11.8548311, 11.8548311, 12.8508802, 12.8508802, 13.1532949, 13.1532949}};

    TestHelper::dumpGolden("correctTemperatureSolution", results.temperature.values);
    TestHelper::expectNear(correctTemperatureSolution, results.temperature.values, 1e-6);

    const std::vector<double> correctTemperatureError{3.60829416e-07, 1.77094312e-07, 1.07265926e-07, 6.91041952e-08, 5.98746289e-08, 4.35682627e-08, 1.60721197e-08, 1.29682699e-08, 4.12447074e-08, 7.01663188e-08};
    TestHelper::dumpGolden("correctTemperatureError", results.temperature.errors);
    TestHelper::expectNear(correctTemperatureError, results.temperature.errors, 1e-6);
}
