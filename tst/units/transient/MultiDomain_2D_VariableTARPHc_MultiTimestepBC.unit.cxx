#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "SlabCreator.hxx"
#include "TestHelpers.hxx"
#include "TestMaterials.hxx"

TEST(MultiDomain_2D_VariableTARPHc_MultiTimestepBC, TestExample_1)
{
    HygroThermFEM::MultiDomain multiDomain;

    const auto & material =
      multiDomain.materials().createSolidMaterial(TestHelper::CottaerSandstone());

    constexpr HygroThermFEM::State state(
      {.temperature = 0.0, .humidity = 0.0, .pressure = 101325.0, .liquidPercent = 1.0});

    TestHelper::SlabBuilder(multiDomain)
      .gridXCoordinates({0, 0.05, 0.1})
      .height(0.05)
      .material(material.name())
      .state(state)
      .startCorner(TestHelper::StartCorner::BottomRight)
      .direction(TestHelper::Direction::CounterClockwise)
      .build();

    /// Create Boundary Conditions

    // Variable boundary conditions (temperature and humidity) over ten timesteps.
    const std::vector<HygroThermFEM::TARPCoefficients> bcCoeff{{20.0, 0.6},
                                                               {20.0, 0.5},
                                                               {20.0, 0.4},
                                                               {20.0, 0.3},
                                                               {20.0, 0.2},
                                                               {18.0, 0.2},
                                                               {16.0, 0.2},
                                                               {14.0, 0.2},
                                                               {12.0, 0.2},
                                                               {10.0, 0.2}};

    const auto surfaceTilt{90.0};

    multiDomain.createBC_TARPHc(1, 2, bcCoeff, surfaceTilt);

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

    std::vector<std::vector<double>> correctWaterContentSolution{{3.46405341, 3.46405341, 0.0033612209, 0.0033612209, 6.52286979e-06, 6.52286979e-06},
 {5.44141426, 5.44141426, 0.0103061278, 0.0103061278, 2.97901105e-05, 2.97901105e-05},
 {6.48364376, 6.48364376, 0.0188496859, 0.0188496859, 7.80855872e-05, 7.80855872e-05},
 {6.77604491, 6.77604491, 0.0282172863, 0.0282172863, 0.000157605642, 0.000157605642},
 {6.45223335, 6.45223335, 0.0377370933, 0.0377370933, 0.000271642682, 0.000271642682},
 {6.02099287, 6.02099287, 0.0471610283, 0.0471610283, 0.000421439041, 0.000421439041},
 {5.51397452, 5.51397452, 0.0563089476, 0.0563089476, 0.000606607468, 0.000606607468},
 {4.99645088, 4.99645088, 0.064819179, 0.064819179, 0.000824628625, 0.000824628625},
 {4.50911743, 4.50911743, 0.0724887629, 0.0724887629, 0.00107157, 0.00107157},
 {4.05915393, 4.05915393, 0.0793325442, 0.0793325442, 0.00134303659, 0.00134303659}};

    TestHelper::dumpGolden("correctWaterContentSolution", waterContentSolution);
    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{{3.96435011, 3.96435011, 2.05915358, 2.05915358, 1.54113049, 1.54113049},
 {5.54103363, 5.54103363, 3.79901598, 3.79901598, 3.23099375, 3.23099375},
 {6.50145244, 6.50145244, 5.12714593, 5.12714593, 4.6501208, 4.6501208},
 {7.15408222, 7.15408222, 6.11626314, 6.11626314, 5.74741213, 5.74741213},
 {7.59524692, 7.59524692, 6.83501096, 6.83501096, 6.56138744, 6.56138744},
 {7.89134988, 7.89134988, 7.34677963, 7.34677963, 7.14918008, 7.14918008},
 {8.0619835, 8.0619835, 7.6913888, 7.6913888, 7.55496526, 7.55496526},
 {8.12822216, 8.12822216, 7.89951153, 7.89951153, 7.81281278, 7.81281278},
 {8.1020303, 8.1020303, 7.99253628, 7.99253628, 7.94730138, 7.94730138},
 {7.98232366, 7.98232366, 7.98057835, 7.98057835, 7.97218498, 7.97218498}};

    TestHelper::dumpGolden("correctTemperatureSolution", temperatureSolution);
    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-6);
}
