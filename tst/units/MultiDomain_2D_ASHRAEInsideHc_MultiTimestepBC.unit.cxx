#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "SlabCreator.hxx"
#include "TestHelpers.hxx"
#include "TestMaterials.hxx"

TEST(MultiDomain_2D_ASHRAEInsideHc_MultiTimestepBC, TestExample_1)
{
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

    // Variable boundary conditions (temperature, pressure and humidity) over ten timesteps.
    const std::vector<HygroThermFEM::ASHRAEInsideCoefficients> bcCoeff{{20.0, 0.6, 101325.0},
                                                                       {20.0, 0.5, 101325.0},
                                                                       {20.0, 0.4, 101325.0},
                                                                       {20.0, 0.3, 101325.0},
                                                                       {20.0, 0.2, 101325.0},
                                                                       {18.0, 0.2, 101325.0},
                                                                       {16.0, 0.2, 101325.0},
                                                                       {14.0, 0.2, 101325.0},
                                                                       {12.0, 0.2, 101325.0},
                                                                       {10.0, 0.2, 101325.0}};

    constexpr auto surfaceTilt{90.0};    // degrees
    constexpr auto surfaceHeight{1.0};   // meters

    multiDomain.createBC_ASHRAEInsideHc(1, 2, bcCoeff, surfaceHeight, surfaceTilt);

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
      {2.914084, 2.914084, 0.003164, 0.003164, 0.000007, 0.000007},
      {4.341478, 4.341478, 0.008334, 0.008334, 0.000027, 0.000027},
      {5.073479, 5.073479, 0.014960, 0.014960, 0.000067, 0.000067},
      {5.399639, 5.399639, 0.022510, 0.022510, 0.000132, 0.000132},
      {5.273585, 5.273585, 0.030437, 0.030437, 0.000228, 0.000228},
      {5.124186, 5.124186, 0.038516, 0.038516, 0.000354, 0.000354},
      {4.928666, 4.928666, 0.046561, 0.046561, 0.000513, 0.000513},
      {4.707886, 4.707886, 0.054427, 0.054427, 0.000704, 0.000704},
      {4.480177, 4.480177, 0.062017, 0.062017, 0.000925, 0.000925},
      {4.265691, 4.265691, 0.069288, 0.069288, 0.001173, 0.001173}};

    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{
      {3.538081, 3.538081, 1.837881, 1.837881, 1.375524, 1.375524},
      {5.112334, 5.112334, 3.477675, 3.477675, 2.948832, 2.948832},
      {6.090571, 6.090571, 4.764796, 4.764796, 4.307946, 4.307946},
      {6.771685, 6.771685, 5.746420, 5.746420, 5.384532, 5.384532},
      {7.248991, 7.248991, 6.478520, 6.478520, 6.203291, 6.203291},
      {7.560374, 7.560374, 7.003452, 7.003452, 6.802140, 6.802140},
      {7.735202, 7.735202, 7.356243, 7.356243, 7.216830, 7.216830},
      {7.797732, 7.797732, 7.566413, 7.566413, 7.478449, 7.478449},
      {7.770254, 7.770254, 7.659923, 7.659923, 7.614250, 7.614250},
      {7.677814, 7.677814, 7.662426, 7.662426, 7.650286, 7.650286}};

    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-6);
}
