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

    std::vector<std::vector<double>> correctWaterContentSolution{
      {3.154369, 3.154369, 0.003407, 0.003407, 0.000007, 0.000007},
      {4.625021, 4.625021, 0.008968, 0.008968, 0.000029, 0.000029},
      {5.384263, 5.384263, 0.016086, 0.016086, 0.000073, 0.000073},
      {5.844827, 5.844827, 0.024160, 0.024160, 0.000145, 0.000145},
      {5.578137, 5.578137, 0.032576, 0.032576, 0.000251, 0.000251},
      {5.231993, 5.231993, 0.041075, 0.041075, 0.000390, 0.000390},
      {4.965266, 4.965266, 0.049438, 0.049438, 0.000565, 0.000565},
      {4.659495, 4.659495, 0.057489, 0.057489, 0.000773, 0.000773},
      {4.329571, 4.329571, 0.065092, 0.065092, 0.001012, 0.001012},
      {3.988121, 3.988121, 0.072146, 0.072146, 0.001280, 0.001280}};

    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{
      {3.995839, 3.995839, 2.075604, 2.075604, 1.553442, 1.553442},
      {5.646619, 5.646619, 3.861483, 3.861483, 3.280844, 3.280844},
      {6.634772, 6.634772, 5.224996, 5.224996, 4.735896, 4.735896},
      {7.306077, 7.306077, 6.240790, 6.240790, 5.862192, 5.862192},
      {7.766054, 7.766054, 6.982390, 6.982390, 6.700567, 6.700567},
      {8.084152, 8.084152, 7.516736, 7.516736, 7.311395, 7.311395},
      {8.276204, 8.276204, 7.883354, 7.883354, 7.739447, 7.739447},
      {8.347525, 8.347525, 8.104696, 8.104696, 8.012790, 8.012790},
      {8.309361, 8.309361, 8.198120, 8.198120, 8.151475, 8.151475},
      {8.171946, 8.171946, 8.177636, 8.177636, 8.171033, 8.171033}};

    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-6);
}
