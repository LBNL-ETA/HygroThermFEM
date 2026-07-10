#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "SlabCreator.hxx"
#include "TestHelpers.hxx"
#include "TestMaterials.hxx"

TEST(MultiDomain_2D_YazdanianKlemsHc_MultiTimestepBC, TestExample_1)
{
    HygroThermFEM::MultiDomain multiDomain;

    const auto & material = multiDomain.materials().createSolidMaterial(TestHelper::CottaerSandstone());

    constexpr HygroThermFEM::State state({
        .temperature = 0.0,
        .humidity = 0.0,
        .pressure = 101325.0,
        .liquidPercent = 1.0
    });

    TestHelper::SlabBuilder(multiDomain)
        .gridXCoordinates({0, 0.05, 0.1})
        .height(0.05)
        .material(material.name())
        .state(state)
        .startCorner(TestHelper::StartCorner::BottomRight)
        .direction(TestHelper::Direction::CounterClockwise)
        .build();

    /// Create Boundary Conditions

    using HygroThermFEM::WindDirection;

    // Variable boundary conditions (temperature, humidity, wind speed and wind direction) over ten timesteps.
    const std::vector<HygroThermFEM::YazdanianKlemsCoefficients> bcCoeff{
      {20.0, 0.6, 3.0, WindDirection::Windward},
      {20.0, 0.5, 3.0, WindDirection::Windward},
      {20.0, 0.4, 3.0, WindDirection::Windward},
      {20.0, 0.3, 4.0, WindDirection::Windward},
      {20.0, 0.2, 4.2, WindDirection::Windward},
      {18.0, 0.2, 4.6, WindDirection::Leeward},
      {16.0, 0.2, 5.0, WindDirection::Leeward},
      {14.0, 0.2, 5.3, WindDirection::Leeward},
      {12.0, 0.2, 5.5, WindDirection::Leeward},
      {10.0, 0.2, 5.9, WindDirection::Leeward}};

    multiDomain.createBC_YazdanianKlemsHc(1, 2, bcCoeff);

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

    std::vector<std::vector<double>> correctWaterContentSolution{{5.82829300489, 5.82829300489, 0.00540581599819, 0.00540581599819, 1.04906624687e-05, 1.04906624687e-05},
 {8.55402847454, 8.55402847454, 0.0161460015962, 0.0161460015962, 5.05647253278e-05, 5.05647253278e-05},
 {9.2638309974, 9.26383099742, 0.029209443204, 0.029209443204, 0.000138165276306, 0.000138165276306},
 {8.47173570575, 8.47173570577, 0.0425423326812, 0.0425423326812, 0.000282968442712, 0.000282968442712},
 {6.5025285266, 6.50252852662, 0.0546797221258, 0.0546797221258, 0.000487944467094, 0.000487944467094},
 {5.0348585616, 5.03485856161, 0.0653204076481, 0.065320407648, 0.000749259884561, 0.000749259884561},
 {4.01370304954, 4.01370304958, 0.0738242381485, 0.0738242381484, 0.00105587655915, 0.00105587655915},
 {3.28355073337, 3.2835507334, 0.0807863518676, 0.0807863518676, 0.00139965090039, 0.00139965090039},
 {2.75847915809, 2.75847915811, 0.0865478712196, 0.0865478712195, 0.00177167465881, 0.00177167465881},
 {2.36802903557, 2.36802903559, 0.0913377665903, 0.0913377665901, 0.00216237664848, 0.00216237664848}};

    TestHelper::dumpGolden("correctWaterContentSolution", waterContentSolution);
    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{{6.46707594886, 6.46707594886, 3.36115344774, 3.36115344774, 2.51558855766, 2.51558855766},
 {8.63726226688, 8.63726226688, 5.99481800129, 5.99481800129, 5.11955279921, 5.11955279921},
 {9.66603030331, 9.66603030331, 7.79244542404, 7.79244542404, 7.12004358095, 7.12004358095},
 {10.6469831774, 10.6469831774, 9.19297992948, 9.19297992948, 8.67152768532, 8.67152768532},
 {11.1727758037, 11.1727758037, 10.1588722332, 10.1588722332, 9.78475988662, 9.78475988662},
 {11.3574179555, 11.3574179555, 10.7377610194, 10.7377610194, 10.4980968216, 10.4980968216},
 {11.5487477435, 11.5487477435, 11.1319501204, 11.1319501204, 10.9725927748, 10.9725927748},
 {11.5334829286, 11.5334829286, 11.3232218701, 11.3232218701, 11.2351317838, 11.2351317838},
 {11.2914200133, 11.2914200133, 11.2981158062, 11.2981158062, 11.282402884, 11.282402884},
 {10.8274795295, 10.8274795295, 11.0540053696, 11.0540053696, 11.1116069688, 11.1116069688}};

    TestHelper::dumpGolden("correctTemperatureSolution", temperatureSolution);
    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-6);
}
