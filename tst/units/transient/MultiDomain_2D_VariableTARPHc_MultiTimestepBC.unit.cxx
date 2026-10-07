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

    std::vector<std::vector<double>> correctWaterContentSolution{{3.46405340748, 3.46405340748, 0.00336122089849, 0.00336122089849, 6.5228697861e-06, 6.5228697861e-06},
 {5.4416941348, 5.4416941348, 0.0103018324263, 0.0103018324263, 2.97811965336e-05, 2.97811965336e-05},
 {6.48451002274, 6.48451002272, 0.01883932731, 0.01883932731, 7.80551744441e-05, 7.80551744441e-05},
 {6.77752915529, 6.77752915529, 0.0282010796398, 0.0282010796398, 0.000157543496687, 0.000157543496687},
 {6.45420782767, 6.45420782766, 0.0377159298122, 0.0377159298122, 0.000271543691866, 0.000271543691866},
 {6.02343031773, 6.02343031772, 0.0471360516013, 0.0471360516013, 0.000421303509185, 0.000421303509185},
 {5.51677886844, 5.51677886842, 0.0562811451572, 0.0562811451572, 0.000606439325259, 0.000606439325259},
 {4.99918175219, 4.99918175217, 0.0647916652557, 0.0647916652557, 0.000824442770311, 0.000824442770311},
 {4.51169713929, 4.51169713929, 0.0724621605573, 0.0724621605574, 0.00107138037422, 0.00107138037422},
 {4.06151866013, 4.06151866012, 0.0793073281661, 0.0793073281661, 0.00134285588587, 0.00134285588587}};

    TestHelper::dumpGolden("correctWaterContentSolution", waterContentSolution);
    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{{3.96182969328, 3.96182969328, 2.06012547793, 2.06012547793, 1.54186093011, 1.54186093011},
 {5.5351839097, 5.5351839097, 3.80106560571, 3.80106560571, 3.23272267104, 3.23272267104},
 {6.49388021835, 6.49388021835, 5.13001798208, 5.13001798208, 4.65272913357, 4.65272913357},
 {7.14576885224, 7.14576885224, 6.11986110497, 6.11986110497, 5.75080132735, 5.75080132735},
 {7.58694268755, 7.58694268755, 6.83919342693, 6.83919342693, 6.56542898143, 6.56542898143},
 {7.88322742513, 7.88322742513, 7.35136823305, 7.35136823305, 7.15370903078, 7.15370903078},
 {8.05437654451, 8.05437654451, 7.69631576141, 7.69631576141, 7.55988926269, 7.55988926269},
 {8.12130654331, 8.12130654331, 7.90454791912, 7.90454791912, 7.81793590724, 7.81793590724},
 {8.09609246131, 8.09609246131, 7.99755521281, 7.99755521281, 7.95247714257, 7.95247714257},
 {7.9773919832, 7.9773919832, 7.98552561458, 7.98552561458, 7.97733336636, 7.97733336636}};

    TestHelper::dumpGolden("correctTemperatureSolution", temperatureSolution);
    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-6);
}
