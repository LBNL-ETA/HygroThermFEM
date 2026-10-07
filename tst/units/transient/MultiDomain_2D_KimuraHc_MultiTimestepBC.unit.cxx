#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "SlabCreator.hxx"
#include "TestHelpers.hxx"
#include "TestMaterials.hxx"

TEST(MultiDomain_2D_KimuraHc_MultiTimestepBC, TestExample_1)
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

    // Variable boundary conditions (temperature, wind speed and humidity) over ten timesteps.
    const std::vector<HygroThermFEM::KimuraCoefficients> bcCoeff{
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

    multiDomain.createBC_KimuraHc(1, 2, bcCoeff);

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

    std::vector<std::vector<double>> correctWaterContentSolution{{8.34404319873, 8.34404319873, 0.00665912024366, 0.00665912024366, 1.2922856212e-05, 1.2922856212e-05},
 {9.57481802244, 9.57481802244, 0.0235678991796, 0.0235678991796, 8.53097121415e-05, 8.53097121415e-05},
 {7.6722403556, 7.6722403556, 0.0430757002183, 0.0430757002183, 0.000272202610423, 0.000272202610423},
 {4.70866261844, 4.70866261844, 0.059997923097, 0.059997923097, 0.000610042603474, 0.000610042603474},
 {2.54283136663, 2.54283136667, 0.0714605467301, 0.07146054673, 0.00112907288867, 0.00112907288867},
 {1.61381538111, 1.61381538114, 0.0805771929858, 0.0805771929856, 0.00188381944486, 0.00188381944486},
 {1.20858382866, 1.20858382868, 0.0873763136349, 0.0873763136347, 0.00282210284749, 0.00282210284749},
 {0.986096707994, 0.986096708001, 0.0928253089951, 0.0928253089949, 0.00390155812244, 0.00390155812244},
 {0.850397739822, 0.850397739825, 0.0972704470086, 0.0972704470084, 0.00507088116443, 0.00507088116443},
 {0.761848597619, 0.761848597621, 0.100918145624, 0.100918145624, 0.00627668226264, 0.00627668226265}};

    TestHelper::dumpGolden("correctWaterContentSolution", waterContentSolution);
    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{{12.2099651321, 12.2099651321, 6.34301202596, 6.34301202596, 4.74729857934, 4.74729857934},
 {16.2922531793, 16.2922531793, 11.3096881565, 11.3096881565, 9.65879509827, 9.65879509827},
 {19.2065553091, 19.2065553091, 15.2076971353, 15.2076971353, 13.8118133325, 13.8118133325},
 {23.1617635712, 23.1617635712, 19.1696685883, 19.1696685883, 17.8219214138, 17.8219214138},
 {27.4161308323, 27.4161308322, 23.2874697517, 23.2874697517, 21.9127345363, 21.9127345363},
 {28.1143743907, 28.1143743907, 25.6235506963, 25.6235506963, 24.6904173697, 24.6904173697},
 {28.6562579234, 28.6562579234, 27.0834763722, 27.0834763722, 26.4819889485, 26.4819889485},
 {28.6080479883, 28.6080479883, 27.8020584036, 27.8020584036, 27.4706124021, 27.4706124021},
 {28.0236535785, 28.0236535785, 27.8780471637, 27.8780471637, 27.7762664638, 27.7762664638},
 {26.9932918941, 26.9932918941, 27.4084645372, 27.4084645372, 27.5017460815, 27.5017460815}};

    TestHelper::dumpGolden("correctTemperatureSolution", temperatureSolution);
    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-6);
}
