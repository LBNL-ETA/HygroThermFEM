#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "SlabCreator.hxx"
#include "TestHelpers.hxx"
#include "TestMaterials.hxx"

TEST(MultiDomain_2D_ASHRAEOutsideHc_MultiTimestepBC, TestExample_1)
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

    // Variable boundary conditions (temperature, humidity and wind speed) over ten timesteps.
    const std::vector<HygroThermFEM::ASHRAEOutsideCoefficients> bcCoeff{{20.0, 0.6, 3},
                                                                       {20.0, 0.5, 3},
                                                                       {20.0, 0.4, 3},
                                                                       {20.0, 0.3, 4},
                                                                       {20.0, 0.2, 4.2},
                                                                       {18.0, 0.2, 4.6},
                                                                       {16.0, 0.2, 5},
                                                                       {14.0, 0.2, 5.3},
                                                                       {12.0, 0.2, 5.5},
                                                                       {10.0, 0.2, 5.9}};

    multiDomain.createBC_ASHRAEOutsideHc(1, 2, bcCoeff);

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
      {8.492737, 8.492737, 0.006610, 0.006610, 0.000013, 0.000013},
      {10.278959, 10.278959, 0.018010, 0.018010, 0.000066, 0.000066},
      {8.868085, 8.868085, 0.031144, 0.031144, 0.000183, 0.000183},
      {6.017858, 6.017858, 0.043230, 0.043230, 0.000372, 0.000372},
      {3.810326, 3.810326, 0.052410, 0.052410, 0.000626, 0.000626},
      {2.814421, 2.814421, 0.059771, 0.059771, 0.000939, 0.000939},
      {2.305688, 2.305688, 0.066020, 0.066020, 0.001304, 0.001304},
      {2.025905, 2.025905, 0.071509, 0.071509, 0.001707, 0.001707},
      {1.856835, 1.856835, 0.076387, 0.076387, 0.002133, 0.002133},
      {1.745478, 1.745478, 0.080724, 0.080724, 0.002565, 0.002565}};

    TestHelper::expectNear(correctWaterContentSolution, waterContentSolution, 1e-6);

    std::vector<std::vector<double>> correctTemperatureSolution{
      {11.860950, 11.860950, 6.037252, 6.037252, 4.422219, 4.422219},
      {12.965216, 12.965216, 9.419692, 9.419692, 8.162457, 8.162457},
      {13.511319, 13.511319, 11.377378, 11.377378, 10.568573, 10.568573},
      {14.478654, 14.478654, 12.880540, 12.880540, 12.298885, 12.298885},
      {15.361766, 15.361766, 14.092252, 14.092252, 13.641057, 13.641057},
      {15.700908, 15.700908, 14.868130, 14.868130, 14.559396, 14.559396},
      {15.380797, 15.380797, 15.093409, 15.093409, 14.959031, 14.959031},
      {14.558736, 14.558736, 14.797427, 14.797427, 14.838052, 14.838052},
      {13.388995, 13.388995, 14.070413, 14.070413, 14.263506, 14.263506},
      {11.942351, 11.942351, 12.989486, 12.989486, 13.309978, 13.309978}};

    TestHelper::expectNear(correctTemperatureSolution, temperatureSolution, 1e-6);
}
