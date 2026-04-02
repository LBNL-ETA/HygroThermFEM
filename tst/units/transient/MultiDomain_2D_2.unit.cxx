#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "SlabCreator.hxx"
#include "TestHelpers.hxx"
#include "TestMaterials.hxx"

using HygroThermFEM::Nodes;

TEST(MultiDomain_2D_2, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Simple two elements example with moisture and heat transfer.");

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
    constexpr auto hc = 1.0;
    constexpr auto airTemperature = 20.0;
    constexpr auto humidity = 0.6;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff{airTemperature, hc, humidity};

    multiDomain.createBC_FixedHc(1, 2, bcCoeff);

    constexpr auto dTime = 3600;
    constexpr auto nSteps = 10;

    auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
    auto humidities = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);
    std::vector<std::vector<double>> temperatureSolution;
    std::vector<std::vector<double>> waterContentSolution;

    for(auto i = 0; i < nSteps; ++i)
    {
        auto aSolution = multiDomain.transient(temperatures, humidities, dTime);
        temperatureSolution.push_back(aSolution.temperature);
        waterContentSolution.push_back(aSolution.waterContent);
        temperatures = aSolution.temperature;
        humidities = aSolution.humidity;
    }

    std::vector<std::vector<double>> correctWaterContentSolution{
      {1.094363, 1.094363, 0.001196, 0.001196, 0.000003, 0.000003},
      {2.123689, 2.123689, 0.003621, 0.003621, 0.000011, 0.000011},
      {3.016116, 3.016116, 0.007207, 0.007207, 0.000028, 0.000028},
      {3.744640, 3.744640, 0.011795, 0.011795, 0.000057, 0.000057},
      {4.361659, 4.361659, 0.017302, 0.017302, 0.000102, 0.000102},
      {4.896016, 4.896016, 0.023670, 0.023670, 0.000166, 0.000166},
      {5.427790, 5.427790, 0.030854, 0.030854, 0.000253, 0.000253},
      {6.240118, 6.240118, 0.038816, 0.038816, 0.000366, 0.000366},
      {6.967339, 6.967339, 0.047525, 0.047525, 0.000510, 0.000510},
      {7.627128, 7.627128, 0.056971, 0.056971, 0.000689, 0.000689}};

    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{
      {1.282432, 1.282432, 0.666356, 0.666356, 0.498720, 0.498720},
      {2.109294, 2.109294, 1.393945, 1.393945, 1.168732, 1.168732},
      {2.808719, 2.808719, 2.099133, 2.099133, 1.865069, 1.865069},
      {3.455716, 3.455716, 2.772728, 2.772728, 2.544384, 2.544384},
      {4.068310, 4.068310, 3.415201, 3.415201, 3.196123, 3.196123},
      {4.651591, 4.651591, 4.027975, 4.027975, 3.818698, 3.818698},
      {5.207521, 5.207521, 4.612350, 4.612350, 4.412679, 4.412679},
      {5.737045, 5.737045, 5.169320, 5.169320, 4.978958, 4.978958},
      {6.241380, 6.241380, 5.700017, 5.700017, 5.518602, 5.518602},
      {6.721448, 6.721448, 6.205448, 6.205448, 6.032637, 6.032637}};

    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-6);
}

TEST(MultiDomain_2D_2, TestExample_1_Repeat)
{
    SCOPED_TRACE("Begin Test: Repeatability test.");

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
    constexpr auto hc = 1.0;
    constexpr auto airTemperature = 20.0;
    constexpr auto humidity = 0.6;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff{airTemperature, hc, humidity};

    multiDomain.createBC_FixedHc(1, 2, bcCoeff);

    constexpr auto dTime = 3600;
    constexpr auto nSteps = 10;

    auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
    auto humidities = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);
    std::vector<std::vector<double>> temperatureSolution;
    std::vector<std::vector<double>> waterContentSolution;

    for(auto i = 0; i < nSteps; ++i)
    {
        auto aSolution = multiDomain.transient(temperatures, humidities, dTime);
        temperatureSolution.push_back(aSolution.temperature);
        waterContentSolution.push_back(aSolution.waterContent);
        temperatures = aSolution.temperature;
        humidities = aSolution.humidity;
    }

    std::vector<std::vector<double>> correctWaterContentSolution{
            {1.094363, 1.094363, 0.001196, 0.001196, 0.000003, 0.000003},
            {2.123689, 2.123689, 0.003621, 0.003621, 0.000011, 0.000011},
            {3.016116, 3.016116, 0.007207, 0.007207, 0.000028, 0.000028},
            {3.744640, 3.744640, 0.011795, 0.011795, 0.000057, 0.000057},
            {4.361659, 4.361659, 0.017302, 0.017302, 0.000102, 0.000102},
            {4.896016, 4.896016, 0.023670, 0.023670, 0.000166, 0.000166},
            {5.427790, 5.427790, 0.030854, 0.030854, 0.000253, 0.000253},
            {6.240118, 6.240118, 0.038816, 0.038816, 0.000366, 0.000366},
            {6.967339, 6.967339, 0.047525, 0.047525, 0.000510, 0.000510},
            {7.627128, 7.627128, 0.056971, 0.056971, 0.000689, 0.000689}};

    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{
            {1.282432, 1.282432, 0.666356, 0.666356, 0.498720, 0.498720},
            {2.109294, 2.109294, 1.393945, 1.393945, 1.168732, 1.168732},
            {2.808719, 2.808719, 2.099133, 2.099133, 1.865069, 1.865069},
            {3.455716, 3.455716, 2.772728, 2.772728, 2.544384, 2.544384},
            {4.068310, 4.068310, 3.415201, 3.415201, 3.196123, 3.196123},
            {4.651591, 4.651591, 4.027975, 4.027975, 3.818698, 3.818698},
            {5.207521, 5.207521, 4.612350, 4.612350, 4.412679, 4.412679},
            {5.737045, 5.737045, 5.169320, 5.169320, 4.978958, 4.978958},
            {6.241380, 6.241380, 5.700017, 5.700017, 5.518602, 5.518602},
            {6.721448, 6.721448, 6.205448, 6.205448, 6.032637, 6.032637}};

    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-6);
}