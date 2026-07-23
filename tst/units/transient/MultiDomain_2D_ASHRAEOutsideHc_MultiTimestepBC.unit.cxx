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
 {12.2930017527, 12.2930017527, 0.0265523868583, 0.0265523868583, 9.46883782938e-05, 9.46883782938e-05},
 {10.3515520896, 10.3515520896, 0.0445100717052, 0.0445100717052, 0.000260421625864, 0.000260421625864},
 {7.09764415272, 7.09764415272, 0.0594217384392, 0.0594217384392, 0.000511865025978, 0.000511865025978},
 {4.17006849691, 4.17006849691, 0.0699416565296, 0.0699416565296, 0.000837595315071, 0.000837595315071},
 {2.96245892768, 2.96245892766, 0.0777689524612, 0.0777689524612, 0.00122944158752, 0.00122944158752},
 {2.37796319667, 2.37796319666, 0.0841500999088, 0.0841500999089, 0.00167688639563, 0.00167688639563},
 {2.06555954241, 2.06555954241, 0.0895204960034, 0.0895204960034, 0.0021611061764, 0.0021611061764},
 {1.88099770524, 1.88099770523, 0.0940976592689, 0.0940976592689, 0.00266073064827, 0.00266073064827},
 {1.76125586291, 1.76125586291, 0.0980143427135, 0.0980143427136, 0.00315650899216, 0.00315650899216}};

    TestHelper::dumpGolden("correctWaterContentSolution", results.moisture.values);
    TestHelper::expectNear(correctWaterContentSolution, results.moisture.values, 1e-6);

    const std::vector<double> correctHumidityError{3.5365268817e-07, 1.51100170388e-08, 1.31362798947e-07, 3.12840593848e-07, 6.02499589492e-07, 5.00637464714e-07, 2.02695216881e-07, 2.17051060396e-09, 1.26386771143e-07, 2.11931307029e-07};
    TestHelper::dumpGolden("correctHumidityError", results.moisture.errors);
    TestHelper::expectNear(correctHumidityError, results.moisture.errors, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{{11.3406394774, 11.3406394774, 5.88952716646, 5.88952716646, 4.40789713764, 4.40789713764},
 {12.3880934172, 12.3880934172, 9.07883793821, 9.07883793821, 7.9037832349, 7.9037832349},
 {12.9475075646, 12.9475075646, 10.9446558718, 10.9446558718, 10.1797178158, 10.1797178158},
 {13.8917711693, 13.8917711693, 12.3850011513, 12.3850011513, 11.8303176801, 11.8303176801},
 {14.8479547732, 14.8479547732, 13.5992194703, 13.5992194703, 13.1543582123, 13.1543582123},
 {15.3535396147, 15.3535396147, 14.4580809545, 14.4580809545, 14.1302851116, 14.1302851116},
 {15.139180218, 15.139180218, 14.7737325326, 14.7737325326, 14.6120804839, 14.6120804839},
 {14.3881998742, 14.3881998742, 14.5562165769, 14.5562165769, 14.5705196668, 14.5705196668},
 {13.2669492274, 13.2669492274, 13.8915579665, 13.8915579665, 14.0626305594, 14.0626305594},
 {11.8559036179, 11.8559036179, 12.8590732164, 12.8590732164, 13.1621229844, 13.1621229844}};

    TestHelper::dumpGolden("correctTemperatureSolution", results.temperature.values);
    TestHelper::expectNear(correctTemperatureSolution, results.temperature.values, 1e-6);

    const std::vector<double> correctTemperatureError{3.59068921633e-07, 1.77551304949e-07, 1.07486806153e-07, 6.91279132497e-08, 5.9793056504e-08, 4.35215510288e-08, 1.60596035796e-08, 1.29652401544e-08, 4.12345432916e-08, 7.01517552439e-08};
    TestHelper::dumpGolden("correctTemperatureError", results.temperature.errors);
    TestHelper::expectNear(correctTemperatureError, results.temperature.errors, 1e-6);
}
