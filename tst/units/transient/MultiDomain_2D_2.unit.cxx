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

    std::vector<std::vector<double>> correctWaterContentSolution{{1.09461711009, 1.09461711009, 0.00106212275082, 0.00106212275082, 2.06118211503e-06, 2.06118211503e-06},
 {2.12454480671, 2.12454480671, 0.00332453223563, 0.00332453223565, 8.83577652521e-06, 8.8357765252e-06},
 {3.08859131346, 3.08859131346, 0.00681626122493, 0.00681626122491, 2.34874781166e-05, 2.34874781166e-05},
 {3.98616922295, 3.98616922285, 0.0115539919642, 0.0115539919643, 4.96329551307e-05, 4.96329551307e-05},
 {4.81710654008, 4.81710653974, 0.0175457192398, 0.0175457192402, 9.13105227987e-05, 9.13105227984e-05},
 {5.59100411526, 5.59100411492, 0.024614798849, 0.0246147988495, 0.000152504121098, 0.000152504121098},
 {6.3232252824, 6.32322528206, 0.0324736339461, 0.0324736339466, 0.000236767714156, 0.000236767714155},
 {7.01348611424, 7.01348611423, 0.0411280947761, 0.0411280947762, 0.000347906575799, 0.000347906575798},
 {7.66182573168, 7.66182573159, 0.05057737343, 0.0505773734299, 0.000489936593705, 0.000489936593704},
 {8.26854058015, 8.26854058021, 0.0608152762943, 0.0608152762941, 0.000667048459842, 0.000667048459841}};

    TestHelper::dumpGolden("correctWaterContentSolution", waterContentSolution);
    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{{1.28164676035, 1.28164676035, 0.666635708606, 0.666635708606, 0.498930765374, 0.498930765374},
 {2.10738142128, 2.10738142128, 1.39454556478, 1.39454556478, 1.16923773433, 1.16923773433},
 {2.80308932481, 2.80308932481, 2.0987644414, 2.0987644414, 1.8649278828, 1.8649278828},
 {3.4409504337, 3.4409504337, 2.76801275607, 2.76801275607, 2.54083206739, 2.54083206739},
 {4.03832544025, 4.03832544026, 3.40145431523, 3.40145431523, 3.18496151692, 3.18496151692},
 {4.60535138128, 4.60535138128, 4.00238669623, 4.00238669623, 3.79676830279, 3.79676830279},
 {5.15050469959, 5.15050469959, 4.57622105406, 4.57622105406, 4.38016436729, 4.38016436729},
 {5.67133244858, 5.67133244858, 5.12424890088, 5.12424890088, 4.93710060517, 4.93710060517},
 {6.1678457675, 6.1678457675, 5.647177959, 5.647177959, 5.46859763049, 5.46859763049},
 {6.64062815225, 6.64062815225, 6.14568765333, 6.14568765333, 5.97542082928, 5.97542082928}};

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

    std::vector<std::vector<double>> correctWaterContentSolution{{1.09461711009, 1.09461711009, 0.00106212275082, 0.00106212275082, 2.06118211503e-06, 2.06118211503e-06},
 {2.12454480671, 2.12454480671, 0.00332453223563, 0.00332453223565, 8.83577652521e-06, 8.8357765252e-06},
 {3.08859131346, 3.08859131346, 0.00681626122493, 0.00681626122491, 2.34874781166e-05, 2.34874781166e-05},
 {3.98616922295, 3.98616922285, 0.0115539919642, 0.0115539919643, 4.96329551307e-05, 4.96329551307e-05},
 {4.81710654008, 4.81710653974, 0.0175457192398, 0.0175457192402, 9.13105227987e-05, 9.13105227984e-05},
 {5.59100411526, 5.59100411492, 0.024614798849, 0.0246147988495, 0.000152504121098, 0.000152504121098},
 {6.3232252824, 6.32322528206, 0.0324736339461, 0.0324736339466, 0.000236767714156, 0.000236767714155},
 {7.01348611424, 7.01348611423, 0.0411280947761, 0.0411280947762, 0.000347906575799, 0.000347906575798},
 {7.66182573168, 7.66182573159, 0.05057737343, 0.0505773734299, 0.000489936593705, 0.000489936593704},
 {8.26854058015, 8.26854058021, 0.0608152762943, 0.0608152762941, 0.000667048459842, 0.000667048459841}};

    TestHelper::dumpGolden("correctWaterContentSolution", waterContentSolution);
    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{{1.28164676035, 1.28164676035, 0.666635708606, 0.666635708606, 0.498930765374, 0.498930765374},
 {2.10738142128, 2.10738142128, 1.39454556478, 1.39454556478, 1.16923773433, 1.16923773433},
 {2.80308932481, 2.80308932481, 2.0987644414, 2.0987644414, 1.8649278828, 1.8649278828},
 {3.4409504337, 3.4409504337, 2.76801275607, 2.76801275607, 2.54083206739, 2.54083206739},
 {4.03832544025, 4.03832544026, 3.40145431523, 3.40145431523, 3.18496151692, 3.18496151692},
 {4.60535138128, 4.60535138128, 4.00238669623, 4.00238669623, 3.79676830279, 3.79676830279},
 {5.15050469959, 5.15050469959, 4.57622105406, 4.57622105406, 4.38016436729, 4.38016436729},
 {5.67133244858, 5.67133244858, 5.12424890088, 5.12424890088, 4.93710060517, 4.93710060517},
 {6.1678457675, 6.1678457675, 5.647177959, 5.647177959, 5.46859763049, 5.46859763049},
 {6.64062815225, 6.64062815225, 6.14568765333, 6.14568765333, 5.97542082928, 5.97542082928}};

    TestHelper::dumpGolden("correctTemperatureSolution", temperatureSolution);
    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-6);
}