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
 {2.12454480847, 2.12454480846, 0.00332453129981, 0.00332453129981, 8.83577459948e-06, 8.83577459948e-06},
 {3.08859132103, 3.08859132102, 0.00681625717886, 0.00681625717887, 2.34874674764e-05, 2.34874674764e-05},
 {3.98616925234, 3.98616925233, 0.0115539760898, 0.0115539760898, 4.96329084918e-05, 4.96329084918e-05},
 {4.81710661998, 4.81710661997, 0.0175456754997, 0.0175456754997, 9.13103720173e-05, 9.13103720173e-05},
 {5.59100058294, 5.59100058293, 0.0246147571232, 0.0246147571233, 0.000152503866564, 0.000152503866564},
 {6.32322206002, 6.32322206, 0.0324734938157, 0.0324734938157, 0.000236767095012, 0.000236767095012},
 {7.01348331165, 7.01348331164, 0.0411277952676, 0.0411277952676, 0.000347905144691, 0.000347905144691},
 {7.66182349256, 7.66182349255, 0.0505768333686, 0.0505768333686, 0.000489933640417, 0.000489933640416},
 {8.26853911198, 8.26853911196, 0.0608143782584, 0.0608143782584, 0.000667042880095, 0.000667042880095}};

    TestHelper::dumpGolden("correctWaterContentSolution", waterContentSolution);
    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{{1.28164676035, 1.28164676035, 0.666635708606, 0.666635708606, 0.498930765374, 0.498930765374},
 {2.10738142121, 2.10738142121, 1.39454556475, 1.39454556475, 1.16923773431, 1.16923773431},
 {2.8030893245, 2.8030893245, 2.09876444125, 2.09876444125, 1.86492788268, 1.86492788268},
 {3.4409504324, 3.4409504324, 2.76801275539, 2.76801275539, 2.54083206683, 2.54083206683},
 {4.03832543633, 4.03832543633, 3.40145431303, 3.40145431303, 3.1849615151, 3.1849615151},
 {4.60535145232, 4.60535145232, 4.00238673122, 4.00238673122, 3.79676832848, 3.79676832848},
 {5.150504815, 5.150504815, 4.57622112893, 4.57622112893, 4.38016442965, 4.38016442965},
 {5.67133259598, 5.67133259598, 5.12424901149, 5.12424901149, 4.93710070334, 4.93710070334},
 {6.16784593627, 6.16784593627, 5.64717809853, 5.64717809853, 5.46859775904, 5.46859775904},
 {6.64062832895, 6.64062832895, 6.14568781204, 6.14568781204, 5.97542097938, 5.97542097938}};

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
 {2.12454480847, 2.12454480846, 0.00332453129981, 0.00332453129981, 8.83577459948e-06, 8.83577459948e-06},
 {3.08859132103, 3.08859132102, 0.00681625717886, 0.00681625717887, 2.34874674764e-05, 2.34874674764e-05},
 {3.98616925234, 3.98616925233, 0.0115539760898, 0.0115539760898, 4.96329084918e-05, 4.96329084918e-05},
 {4.81710661998, 4.81710661997, 0.0175456754997, 0.0175456754997, 9.13103720173e-05, 9.13103720173e-05},
 {5.59100058294, 5.59100058293, 0.0246147571232, 0.0246147571233, 0.000152503866564, 0.000152503866564},
 {6.32322206002, 6.32322206, 0.0324734938157, 0.0324734938157, 0.000236767095012, 0.000236767095012},
 {7.01348331165, 7.01348331164, 0.0411277952676, 0.0411277952676, 0.000347905144691, 0.000347905144691},
 {7.66182349256, 7.66182349255, 0.0505768333686, 0.0505768333686, 0.000489933640417, 0.000489933640416},
 {8.26853911198, 8.26853911196, 0.0608143782584, 0.0608143782584, 0.000667042880095, 0.000667042880095}};

    TestHelper::dumpGolden("correctWaterContentSolution", waterContentSolution);
    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{{1.28164676035, 1.28164676035, 0.666635708606, 0.666635708606, 0.498930765374, 0.498930765374},
 {2.10738142121, 2.10738142121, 1.39454556475, 1.39454556475, 1.16923773431, 1.16923773431},
 {2.8030893245, 2.8030893245, 2.09876444125, 2.09876444125, 1.86492788268, 1.86492788268},
 {3.4409504324, 3.4409504324, 2.76801275539, 2.76801275539, 2.54083206683, 2.54083206683},
 {4.03832543633, 4.03832543633, 3.40145431303, 3.40145431303, 3.1849615151, 3.1849615151},
 {4.60535145232, 4.60535145232, 4.00238673122, 4.00238673122, 3.79676832848, 3.79676832848},
 {5.150504815, 5.150504815, 4.57622112893, 4.57622112893, 4.38016442965, 4.38016442965},
 {5.67133259598, 5.67133259598, 5.12424901149, 5.12424901149, 4.93710070334, 4.93710070334},
 {6.16784593627, 6.16784593627, 5.64717809853, 5.64717809853, 5.46859775904, 5.46859775904},
 {6.64062832895, 6.64062832895, 6.14568781204, 6.14568781204, 5.97542097938, 5.97542097938}};

    TestHelper::dumpGolden("correctTemperatureSolution", temperatureSolution);
    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-6);
}