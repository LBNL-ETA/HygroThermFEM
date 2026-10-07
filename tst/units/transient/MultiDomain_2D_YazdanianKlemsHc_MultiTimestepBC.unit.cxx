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
 {8.55426817487, 8.55426817487, 0.0161327475916, 0.0161327475916, 5.05280431579e-05, 5.05280431579e-05},
 {9.26471106287, 9.26471106287, 0.0291833416826, 0.0291833416826, 0.000138038344806, 0.000138038344806},
 {8.47346219055, 8.47346219055, 0.0425084940561, 0.0425084940561, 0.000282709662037, 0.000282709662037},
 {6.50467375496, 6.50467375496, 0.0546417773645, 0.0546417773645, 0.00048752698655, 0.00048752698655},
 {5.03659989671, 5.03659989671, 0.0652822217664, 0.0652822217664, 0.000748676656134, 0.000748676656134},
 {4.01508072587, 4.01508072587, 0.0737866318258, 0.0737866318258, 0.00105512584641, 0.00105512584641},
 {3.28458836694, 3.28458836692, 0.0807495301101, 0.0807495301101, 0.00139873377229, 0.00139873377228},
 {2.75923645373, 2.75923645371, 0.0865118115952, 0.0865118115952, 0.001770595373, 0.001770595373},
 {2.36857115253, 2.36857115251, 0.0913023468747, 0.0913023468748, 0.00216114200973, 0.00216114200973}};

    TestHelper::dumpGolden("correctWaterContentSolution", waterContentSolution);
    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{{6.46597006261, 6.46597006261, 3.36118422479, 3.36118422479, 2.51561203167, 2.51561203167},
 {8.63423771626, 8.63423771626, 5.99429108377, 5.99429108377, 5.11916684859, 5.11916684859},
 {9.66194949903, 9.66194949903, 7.79161354216, 7.79161354216, 7.11933104259, 7.11933104259},
 {10.6427738089, 10.6427738089, 9.19229666589, 9.19229666589, 8.67085136091, 8.67085136091},
 {11.1688309071, 11.1688309071, 10.1583742216, 10.1583742216, 9.78423960504, 9.78423960504},
 {11.3536183671, 11.3536183671, 10.7372871941, 10.7372871941, 10.4976419894, 10.4976419894},
 {11.5455460288, 11.5455460288, 11.1315324601, 11.1315324601, 10.9722033866, 10.9722033866},
 {11.5308737976, 11.5308737976, 11.3228964777, 11.3228964777, 11.2348333637, 11.2348333637},
 {11.2893569601, 11.2893569601, 11.2979077849, 11.2979077849, 11.2822188322, 11.2822188322},
 {10.825922623, 10.825922623, 11.0539318148, 11.0539318148, 11.111553916, 11.111553916}};

    TestHelper::dumpGolden("correctTemperatureSolution", temperatureSolution);
    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-6);
}
