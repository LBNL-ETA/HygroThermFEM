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
 {12.2929215508, 12.2929215508, 0.026611347293, 0.026611347293, 9.48860366919e-05, 9.48860366919e-05},
 {10.3514766669, 10.3514766669, 0.0445907632394, 0.0445907632394, 0.000260947923744, 0.000260947923744},
 {7.09759401081, 7.09759401081, 0.0595089244201, 0.0595089244201, 0.000512777756459, 0.000512777756459},
 {4.17004866532, 4.17004866532, 0.0700311580521, 0.0700311580521, 0.00083893608344, 0.00083893608344},
 {2.96245003568, 2.96245003568, 0.0778593993037, 0.0778593993037, 0.00123124594344, 0.00123124594344},
 {2.3779590486, 2.3779590486, 0.0842406559254, 0.0842406559254, 0.00167917587333, 0.00167917587333},
 {2.06555751984, 2.06555751984, 0.0896106441001, 0.0896106441001, 0.00216388194817, 0.00216388194817},
 {1.88099661718, 1.88099661718, 0.0941873419864, 0.0941873419864, 0.00266397910339, 0.00266397910339},
 {1.76125514073, 1.76125514073, 0.0981037650848, 0.0981037650848, 0.00316020746145, 0.00316020746145}};

    TestHelper::dumpGolden("correctWaterContentSolution", results.moisture.values);
    TestHelper::expectNear(correctWaterContentSolution, results.moisture.values, 1e-6);

    const std::vector<double> correctHumidityError{3.53642790423e-07, 1.51114339004e-08, 1.31364165813e-07, 3.12840133226e-07, 6.02500729581e-07, 5.0063240943e-07, 2.02691482664e-07, 2.17237657918e-09, 1.26387051926e-07, 2.11930447728e-07};
    TestHelper::dumpGolden("correctHumidityError", results.moisture.errors);
    TestHelper::expectNear(correctHumidityError, results.moisture.errors, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{{11.3406394774, 11.3406394774, 5.88952716646, 5.88952716646, 4.40789713764, 4.40789713764},
 {12.3881074985, 12.3881074985, 9.0788448678, 9.0788448678, 7.90378845265, 7.90378845265},
 {12.947531981, 12.947531981, 10.944671254, 10.944671254, 10.1797307395, 10.1797307395},
 {13.8917957321, 13.8917957321, 12.3850205922, 12.3850205922, 11.8303356273, 11.8303356273},
 {14.8479793289, 14.8479793289, 13.5992409842, 13.5992409842, 13.1543790091, 13.1543790091},
 {15.3535568944, 15.3535568944, 14.4580998293, 14.4580998293, 14.1303046803, 14.1303046803},
 {15.1391921296, 15.1391921296, 14.7737476076, 14.7737476076, 14.612096925, 14.612096925},
 {14.3882080951, 14.3882080951, 14.5562280774, 14.5562280774, 14.5705326601, 14.5705326601},
 {13.2669550519, 13.2669550519, 13.8915665975, 13.8915665975, 14.0626405403, 14.0626405403},
 {11.8559078652, 11.8559078652, 12.8590797013, 12.8590797013, 13.1621305948, 13.1621305948}};

    TestHelper::dumpGolden("correctTemperatureSolution", results.temperature.values);
    TestHelper::expectNear(correctTemperatureSolution, results.temperature.values, 1e-6);

    const std::vector<double> correctTemperatureError{3.59068921633e-07, 1.77551386549e-07, 1.07487040418e-07, 6.91281111031e-08, 5.97931541559e-08, 4.3521390913e-08, 1.60593915461e-08, 1.29654313791e-08, 4.12346949592e-08, 7.01518693871e-08};
    TestHelper::dumpGolden("correctTemperatureError", results.temperature.errors);
    TestHelper::expectNear(correctTemperatureError, results.temperature.errors, 1e-6);
}
