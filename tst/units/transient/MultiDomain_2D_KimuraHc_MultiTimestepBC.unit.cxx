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
 {9.57473720364, 9.57473720364, 0.0236292151785, 0.0236292151785, 8.55215169494e-05, 8.55215169494e-05},
 {7.6721331986, 7.6721331986, 0.0431861207474, 0.0431861207474, 0.000272953789135, 0.000272953789135},
 {4.70860360684, 4.70860360684, 0.0601348562618, 0.0601348562618, 0.000611638717298, 0.000611638717298},
 {2.54279992009, 2.54279992014, 0.0716140169474, 0.0716140169473, 0.00113188355728, 0.00113188355728},
 {1.61379447563, 1.61379447565, 0.080743435298, 0.0807434352979, 0.00188832869238, 0.00188832869238},
 {1.20857344424, 1.20857344426, 0.0875447839843, 0.0875447839842, 0.00282849525249, 0.00282849525249},
 {0.98609178775, 0.986091787758, 0.0929931349256, 0.0929931349254, 0.00390992819263, 0.00390992819263},
 {0.85039545128, 0.850395451284, 0.0974365299717, 0.0974365299715, 0.00508124217706, 0.00508124217706},
 {0.761847530482, 0.761847530484, 0.101082224466, 0.101082224466, 0.00628898114167, 0.00628898114167}};

    TestHelper::dumpGolden("correctWaterContentSolution", waterContentSolution);
    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{{12.2099651321, 12.2099651321, 6.34301202596, 6.34301202596, 4.74729857934, 4.74729857934},
 {16.29227118, 16.29227118, 11.3096969189, 11.3096969189, 9.65880167625, 9.65880167625},
 {19.2065951856, 19.2065951856, 15.2077207865, 15.2077207865, 13.8118328261, 13.8118328261},
 {23.1618257419, 23.1618257419, 19.1697103179, 19.1697103179, 17.8219578421, 17.8219578421},
 {27.4161895564, 27.4161895564, 23.2875179048, 23.2875179049, 21.912780238, 21.912780238},
 {28.1144312037, 28.1144312036, 25.6236014471, 25.6236014471, 24.6904676655, 24.6904676655},
 {28.6563062837, 28.6563062837, 27.0835243599, 27.0835243599, 26.4820385417, 26.4820385417},
 {28.6080880727, 28.6080880727, 27.8021011727, 27.8021011727, 27.4706580436, 27.4706580436},
 {28.0236867909, 28.0236867909, 27.8780841823, 27.8780841823, 27.7763068657, 27.7763068657},
 {26.9933195948, 26.9933195948, 27.4084961743, 27.4084961743, 27.5017811351, 27.5017811351}};

    TestHelper::dumpGolden("correctTemperatureSolution", temperatureSolution);
    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-6);
}
