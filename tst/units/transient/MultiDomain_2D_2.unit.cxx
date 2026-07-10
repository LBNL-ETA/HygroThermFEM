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

    std::vector<std::vector<double>> correctWaterContentSolution{{1.09461711, 1.09461711, 0.00106212274, 0.00106212274, 2.0611821e-06, 2.0611821e-06},
 {2.12453856, 2.12453856, 0.00332476594, 0.00332476594, 8.83614102e-06, 8.83614102e-06},
 {3.08856273, 3.08856273, 0.00681715717, 0.00681715717, 2.34892206e-05, 2.34892206e-05},
 {3.98609241, 3.98609241, 0.0115561183, 0.0115561183, 4.96378973e-05, 4.96378973e-05},
 {4.81694645, 4.81694645, 0.0175497905, 0.0175497905, 9.13214799e-05, 9.13214799e-05},
 {5.59070789, 5.59070789, 0.0246218215, 0.0246218215, 0.00015252549, 0.00015252549},
 {6.32275606, 6.32275606, 0.0324843085, 0.0324843085, 0.000236804918, 0.000236804918},
 {7.01280645, 7.01280645, 0.0411431025, 0.0411431025, 0.000347966351, 0.000347966351},
 {7.66089732, 7.66089732, 0.0505973932, 0.0505973932, 0.000490026515, 0.000490026515},
 {8.26732389, 8.26732389, 0.0608409759, 0.0608409759, 0.000667176594, 0.000667176594}};

    TestHelper::dumpGolden("correctWaterContentSolution", waterContentSolution);
    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{{1.28240792, 1.28240792, 0.666358116, 0.666358116, 0.498722082, 0.498722082},
 {2.10922744, 2.10922744, 1.39394071, 1.39394071, 1.16872946, 1.16872946},
 {2.80614433, 2.80614433, 2.09781943, 2.09781943, 1.86408595, 1.86408595},
 {3.44527441, 3.44527441, 2.7667124, 2.7667124, 2.53963473, 2.53963473},
 {4.04394943, 4.04394943, 3.39977812, 3.39977812, 3.18338587, 3.18338587},
 {4.6121474, 4.6121474, 4.00035564, 4.00035564, 3.79482214, 3.79482214},
 {5.15822554, 5.15822554, 4.57389429, 4.57389429, 4.37789194, 4.37789194},
 {5.67992602, 5.67992602, 5.12161198, 5.12161198, 4.93450008, 4.93450008},
 {6.17729243, 6.17729243, 5.64420648, 5.64420648, 5.46564798, 5.46564798},
 {6.65090903, 6.65090903, 6.14235613, 6.14235613, 5.97209526, 5.97209526}};

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

    std::vector<std::vector<double>> correctWaterContentSolution{{1.09461711, 1.09461711, 0.00106212274, 0.00106212274, 2.0611821e-06, 2.0611821e-06},
 {2.12453856, 2.12453856, 0.00332476594, 0.00332476594, 8.83614102e-06, 8.83614102e-06},
 {3.08856273, 3.08856273, 0.00681715717, 0.00681715717, 2.34892206e-05, 2.34892206e-05},
 {3.98609241, 3.98609241, 0.0115561183, 0.0115561183, 4.96378973e-05, 4.96378973e-05},
 {4.81694645, 4.81694645, 0.0175497905, 0.0175497905, 9.13214799e-05, 9.13214799e-05},
 {5.59070789, 5.59070789, 0.0246218215, 0.0246218215, 0.00015252549, 0.00015252549},
 {6.32275606, 6.32275606, 0.0324843085, 0.0324843085, 0.000236804918, 0.000236804918},
 {7.01280645, 7.01280645, 0.0411431025, 0.0411431025, 0.000347966351, 0.000347966351},
 {7.66089732, 7.66089732, 0.0505973932, 0.0505973932, 0.000490026515, 0.000490026515},
 {8.26732389, 8.26732389, 0.0608409759, 0.0608409759, 0.000667176594, 0.000667176594}};

    TestHelper::dumpGolden("correctWaterContentSolution", waterContentSolution);
    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{{1.28240792, 1.28240792, 0.666358116, 0.666358116, 0.498722082, 0.498722082},
 {2.10922744, 2.10922744, 1.39394071, 1.39394071, 1.16872946, 1.16872946},
 {2.80614433, 2.80614433, 2.09781943, 2.09781943, 1.86408595, 1.86408595},
 {3.44527441, 3.44527441, 2.7667124, 2.7667124, 2.53963473, 2.53963473},
 {4.04394943, 4.04394943, 3.39977812, 3.39977812, 3.18338587, 3.18338587},
 {4.6121474, 4.6121474, 4.00035564, 4.00035564, 3.79482214, 3.79482214},
 {5.15822554, 5.15822554, 4.57389429, 4.57389429, 4.37789194, 4.37789194},
 {5.67992602, 5.67992602, 5.12161198, 5.12161198, 4.93450008, 4.93450008},
 {6.17729243, 6.17729243, 5.64420648, 5.64420648, 5.46564798, 5.46564798},
 {6.65090903, 6.65090903, 6.14235613, 6.14235613, 5.97209526, 5.97209526}};

    TestHelper::dumpGolden("correctTemperatureSolution", temperatureSolution);
    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-6);
}