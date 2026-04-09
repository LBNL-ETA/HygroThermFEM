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

    std::vector<std::vector<double>> correctWaterContentSolution{
      {8.184325, 8.184325, 0.006793, 0.006793, 0.000015, 0.000015},
      {10.105295, 10.105295, 0.018206, 0.018206, 0.000069, 0.000069},
      {8.725453, 8.725453, 0.031337, 0.031337, 0.000188, 0.000188},
      {5.924715, 5.924715, 0.043412, 0.043412, 0.000378, 0.000378},
      {3.780085, 3.780085, 0.052575, 0.052575, 0.000635, 0.000635},
      {2.797976, 2.797976, 0.059927, 0.059927, 0.000951, 0.000951},
      {2.296897, 2.296897, 0.066174, 0.066174, 0.001318, 0.001318},
      {2.020731, 2.020731, 0.071663, 0.071663, 0.001723, 0.001723},
      {1.853538, 1.853538, 0.076541, 0.076541, 0.002151, 0.002151},
      {1.743266, 1.743266, 0.080878, 0.080878, 0.002585, 0.002585}};

    TestHelper::expectNear(correctWaterContentSolution, results.moisture.values, 1e-6);

    const std::vector<double> correctHumidityError{4.940245e-07, 4.597750e-08, 1.193402e-07, 3.399188e-07, 6.045441e-07, 4.100176e-07, 1.502028e-07, 3.533361e-08, 1.475798e-07, 2.267924e-07};
    TestHelper::expectNear(correctHumidityError, results.moisture.errors, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{
      {11.872402, 11.872402, 6.162407, 6.162407, 4.612121, 4.612121},
      {13.069614, 13.069614, 9.542635, 9.542635, 8.302245, 8.302245},
      {13.600559, 13.600559, 11.485074, 11.485074, 10.684342, 10.684342},
      {14.553572, 14.553572, 12.972296, 12.972296, 12.396682, 12.396682},
      {15.424064, 15.424064, 14.169506, 14.169506, 13.723479, 13.723479},
      {15.747747, 15.747747, 14.930261, 14.930261, 14.626632, 14.626632},
      {15.414469, 15.414469, 15.141423, 15.141423, 15.011881, 15.011881},
      {14.582974, 14.582974, 14.833720, 14.833720, 14.878510, 14.878510},
      {13.406563, 13.406563, 14.097522, 14.097522, 14.293973, 14.293973},
      {11.954915, 11.954915, 13.009476, 13.009476, 13.332604, 13.332604}};

    TestHelper::expectNear(correctTemperatureSolution, results.temperature.values, 1e-6);

    const std::vector<double> correctTemperatureError{4.605507e-07, 1.811901e-07, 1.083029e-07, 7.053318e-08, 5.837855e-08, 3.781820e-08, 1.038893e-08, 1.770228e-08, 4.496964e-08, 7.310389e-08};
    TestHelper::expectNear(correctTemperatureError, results.temperature.errors, 1e-6);
}
