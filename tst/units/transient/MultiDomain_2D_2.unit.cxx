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

    std::vector<std::vector<double>> correctWaterContentSolution{{1.09461711009, 1.09461711009, 0.00106212274278, 0.00106212274278, 2.06118209937e-06, 2.06118209937e-06},
 {2.12454458186, 2.12454458185, 0.00332465075882, 0.00332465075883, 8.83604535997e-06, 8.83604535997e-06},
 {3.08859064578, 3.08859064577, 0.00681662029248, 0.00681662029249, 2.34886239931e-05, 2.34886239931e-05},
 {3.98616801359, 3.98616801356, 0.0115546579773, 0.0115546579774, 4.96358079088e-05, 4.96358079087e-05},
 {4.81710476838, 4.81710476835, 0.0175467212214, 0.0175467212214, 9.13160526607e-05, 9.13160526606e-05},
 {5.59099804797, 5.59099804794, 0.0246161844676, 0.0246161844677, 0.000152513487471, 0.000152513487471},
 {6.32321887721, 6.32321887718, 0.032475297321, 0.0324752973211, 0.000236781876674, 0.000236781876674},
 {7.01347951402, 7.01347951397, 0.0411299704612, 0.0411299704612, 0.000347926374135, 0.000347926374135},
 {7.6618191185, 7.66181911845, 0.0505793728134, 0.0505793728134, 0.000489962656491, 0.000489962656491},
 {8.26853420568, 8.26853420565, 0.0608172711673, 0.0608172711674, 0.000667081055937, 0.000667081055937}};

    TestHelper::dumpGolden("correctWaterContentSolution", waterContentSolution);
    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{{1.28164676035, 1.28164676035, 0.666635708606, 0.666635708606, 0.498930765374, 0.498930765374},
 {2.10738142872, 2.10738142872, 1.39454556833, 1.39454556833, 1.16923773706, 1.16923773706},
 {2.80308935283, 2.80308935283, 2.09876445655, 2.09876445655, 1.86492789506, 1.86492789506},
 {3.44095049666, 3.44095049666, 2.76801279377, 2.76801279377, 2.54083209917, 2.54083209917},
 {4.03832555204, 4.03832555204, 3.40145438772, 3.40145438772, 3.18496157997, 3.18496157997},
 {4.60535158975, 4.60535158975, 4.00238683381, 4.00238683381, 3.79676842282, 3.79676842282},
 {5.15050499226, 5.15050499226, 4.57622126484, 4.57622126484, 4.38016455678, 4.38016455678},
 {5.67133282292, 5.67133282292, 5.12424918815, 5.12424918815, 4.93710086971, 4.93710086971},
 {6.16784622061, 6.16784622061, 5.64717832337, 5.64717832337, 5.46859797187, 5.46859797187},
 {6.6406286774, 6.6406286774, 6.14568809211, 6.14568809211, 5.97542124579, 5.97542124579}};

    TestHelper::dumpGolden("correctTemperatureSolution", temperatureSolution);
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

    std::vector<std::vector<double>> correctWaterContentSolution{{1.09461711009, 1.09461711009, 0.00106212274278, 0.00106212274278, 2.06118209937e-06, 2.06118209937e-06},
 {2.12454458186, 2.12454458185, 0.00332465075882, 0.00332465075883, 8.83604535997e-06, 8.83604535997e-06},
 {3.08859064578, 3.08859064577, 0.00681662029248, 0.00681662029249, 2.34886239931e-05, 2.34886239931e-05},
 {3.98616801359, 3.98616801356, 0.0115546579773, 0.0115546579774, 4.96358079088e-05, 4.96358079087e-05},
 {4.81710476838, 4.81710476835, 0.0175467212214, 0.0175467212214, 9.13160526607e-05, 9.13160526606e-05},
 {5.59099804797, 5.59099804794, 0.0246161844676, 0.0246161844677, 0.000152513487471, 0.000152513487471},
 {6.32321887721, 6.32321887718, 0.032475297321, 0.0324752973211, 0.000236781876674, 0.000236781876674},
 {7.01347951402, 7.01347951397, 0.0411299704612, 0.0411299704612, 0.000347926374135, 0.000347926374135},
 {7.6618191185, 7.66181911845, 0.0505793728134, 0.0505793728134, 0.000489962656491, 0.000489962656491},
 {8.26853420568, 8.26853420565, 0.0608172711673, 0.0608172711674, 0.000667081055937, 0.000667081055937}};

    TestHelper::dumpGolden("correctWaterContentSolution", waterContentSolution);
    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{{1.28164676035, 1.28164676035, 0.666635708606, 0.666635708606, 0.498930765374, 0.498930765374},
 {2.10738142872, 2.10738142872, 1.39454556833, 1.39454556833, 1.16923773706, 1.16923773706},
 {2.80308935283, 2.80308935283, 2.09876445655, 2.09876445655, 1.86492789506, 1.86492789506},
 {3.44095049666, 3.44095049666, 2.76801279377, 2.76801279377, 2.54083209917, 2.54083209917},
 {4.03832555204, 4.03832555204, 3.40145438772, 3.40145438772, 3.18496157997, 3.18496157997},
 {4.60535158975, 4.60535158975, 4.00238683381, 4.00238683381, 3.79676842282, 3.79676842282},
 {5.15050499226, 5.15050499226, 4.57622126484, 4.57622126484, 4.38016455678, 4.38016455678},
 {5.67133282292, 5.67133282292, 5.12424918815, 5.12424918815, 4.93710086971, 4.93710086971},
 {6.16784622061, 6.16784622061, 5.64717832337, 5.64717832337, 5.46859797187, 5.46859797187},
 {6.6406286774, 6.6406286774, 6.14568809211, 6.14568809211, 5.97542124579, 5.97542124579}};

    TestHelper::dumpGolden("correctTemperatureSolution", temperatureSolution);
    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-6);
}