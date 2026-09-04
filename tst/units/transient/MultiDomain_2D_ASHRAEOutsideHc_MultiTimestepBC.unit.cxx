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

    std::vector<std::vector<double>> correctWaterContentSolution{{11.6503225665, 11.6503225665, 0.00809627880971, 0.00809627880971, 1.57118422778e-05, 1.57118422778e-05},
 {12.2930012318, 12.2930012318, 0.0265523097984, 0.0265523097985, 9.46881481948e-05, 9.46881481947e-05},
 {10.3515516382, 10.3515516382, 0.0445095707492, 0.0445095707493, 0.000260419529335, 0.000260419529335},
 {7.0976436536, 7.0976436536, 0.0594212378878, 0.0594212378878, 0.000511860813953, 0.000511860813953},
 {4.17006794545, 4.17006794545, 0.0699411555234, 0.069941155523, 0.000837588777057, 0.000837588777057},
 {2.96245873568, 2.96245873596, 0.0777684057663, 0.0777684057678, 0.00122943230916, 0.00122943230915},
 {2.37796311536, 2.37796311554, 0.08414954952, 0.0841495495203, 0.00167687420876, 0.00167687420876},
 {2.0655595014, 2.06555950148, 0.089519947313, 0.0895199473133, 0.00216109104212, 0.00216109104212},
 {1.88099768123, 1.88099768128, 0.094097112963, 0.0940971129629, 0.00266071263541, 0.00266071263541},
 {1.76125584738, 1.7612558474, 0.0980137989028, 0.0980137989029, 0.00315648825176, 0.00315648825176}};

    TestHelper::dumpGolden("correctWaterContentSolution", results.moisture.values);
    TestHelper::expectNear(correctWaterContentSolution, results.moisture.values, 1e-6);

    const std::vector<double> correctHumidityError{3.53652608549e-07, 1.51098026647e-08, 1.31362707695e-07, 3.1284040408e-07, 6.02498628241e-07, 5.00637369386e-07, 2.02695163378e-07, 2.17054146241e-09, 1.26386792257e-07, 2.11931325616e-07};
    TestHelper::dumpGolden("correctHumidityError", results.moisture.errors);
    TestHelper::expectNear(correctHumidityError, results.moisture.errors, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{{11.3406398319, 11.3406398319, 5.88952735183, 5.88952735183, 4.40789727638, 4.40789727638},
 {12.3880935515, 12.3880935515, 9.07883809169, 9.07883809169, 7.90378338463, 7.90378338463},
 {12.9475077527, 12.9475077527, 10.944656045, 10.944656045, 10.1797179825, 10.1797179825},
 {13.8917714277, 13.8917714277, 12.3850013698, 12.3850013698, 11.8303178847, 11.8303178847},
 {14.8479553433, 14.8479553433, 13.5992198709, 13.5992198709, 13.1543585625, 13.1543585625},
 {15.3535399662, 15.353539966, 14.4580813248, 14.4580813248, 14.1302854756, 14.1302854756},
 {15.1391804494, 15.1391804492, 14.7737328314, 14.7737328315, 14.6120807977, 14.6120807977},
 {14.3882000337, 14.3882000336, 14.5562168064, 14.5562168064, 14.570519916, 14.570519916},
 {13.2669493405, 13.2669493405, 13.8915581388, 13.8915581388, 14.0626307495, 14.0626307495},
 {11.8559036978, 11.8559036978, 12.8590733432, 12.8590733432, 13.1621231257, 13.1621231257}};

    TestHelper::dumpGolden("correctTemperatureSolution", results.temperature.values);
    TestHelper::expectNear(correctTemperatureSolution, results.temperature.values, 1e-6);

    const std::vector<double> correctTemperatureError{3.59068986954e-07, 1.77551304898e-07, 1.07486805512e-07, 6.91279142891e-08, 5.97930646169e-08, 4.3521548768e-08, 1.60595995627e-08, 1.29652435638e-08, 4.12345461906e-08, 7.01517579223e-08};
    TestHelper::dumpGolden("correctTemperatureError", results.temperature.errors);
    TestHelper::expectNear(correctTemperatureError, results.temperature.errors, 1e-6);
}
