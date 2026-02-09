#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "SlabCreator.hxx"
#include "TestMaterials.hxx"

using HygroThermFEM::Nodes;
using HygroThermFEM::SimulationProperties;

class MultiDomain_2D_ExcludeLiquidTransport_1 : public testing::Test
{
protected:
    void SetUp() override
    {
        const auto relaxationParameter{0.8};
        const auto errorTolerance{1e-5};
        const auto numberOfIterations{20u};
        SimulationProperties::Instance().setIterationParameters(
          relaxationParameter, errorTolerance, numberOfIterations);
    }

    void TearDown() override
    {
        SimulationProperties::Instance().reset();
    }
};

TEST_F(MultiDomain_2D_ExcludeLiquidTransport_1, TestExample_1)
{
    const auto excludeWaterLiquidTransportation{true};
    const auto excludeHeatOfEvaporation{false};
    const auto excludeCapillaryConduction{false};
    const auto excludeVaporDiffusionConduction{false};
    const auto thermalConductivityMoistureAndTemperatureDependent{false};

    SimulationProperties::Instance().setCalculationParameters(
      excludeWaterLiquidTransportation,
      excludeHeatOfEvaporation,
      excludeCapillaryConduction,
      excludeVaporDiffusionConduction,
      thermalConductivityMoistureAndTemperatureDependent);

    HygroThermFEM::MultiDomain multiDomain;

    const auto & material = multiDomain.materials().createSolidMaterial(TestHelper::CottaerSandstone());

    constexpr HygroThermFEM::State state({
        .temperature = 0.0,
        .humidity = 0.99,
        .pressure = 101325.0,
        .liquidPercent = 1.0
    });

    TestHelper::SlabBuilder(multiDomain)
        .gridXCoordinates({0, 0.05, 0.1})
        .height(0.05)
        .material(material.name())
        .state(state)
        .startCorner(TestHelper::StartCorner::BottomRight)
        .direction(TestHelper::Direction::CounterClockwise)
        .build();

    /// Create Boundary Conditions
    constexpr auto hc = 5.0;
    constexpr auto airTemperature = 10.0;
    constexpr auto humidity = 0.0;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff{airTemperature, hc, humidity};

    multiDomain.createBC_FixedHc(1, 2, bcCoeff);
    multiDomain.createBC_FixedHc(5, 6, bcCoeff);

    constexpr auto dTime = 3600;
    constexpr auto nSteps = 24;

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
      {62.28723, 62.28723, 63.000009, 63.000009, 62.28723, 62.28723},
      {61.553636, 61.553636, 63.0, 63.0, 61.553636, 61.553636},
      {60.798374, 60.798374, 62.999991, 62.999991, 60.798374, 60.798374},
      {60.020186, 60.020186, 62.999976, 62.999976, 60.020186, 60.020186},
      {59.21748, 59.21748, 62.999955, 62.999955, 59.21748, 59.21748},
      {58.388297, 58.388297, 62.999927, 62.999927, 58.388297, 58.388297},
      {57.53025, 57.53025, 62.999893, 62.999893, 57.53025, 57.53025},
      {56.64044, 56.64044, 62.999853, 62.999853, 56.64044, 56.64044},
      {55.71532, 55.71532, 62.999805, 62.999805, 55.71532, 55.71532},
      {54.750515, 54.750515, 62.99975, 62.99975, 54.750515, 54.750515},
      {53.740556, 53.740556, 62.999687, 62.999687, 53.740556, 53.740556},
      {52.678483, 52.678483, 62.999616, 62.999616, 52.678483, 52.678483},
      {51.555228, 51.555228, 62.999537, 62.999537, 51.555228, 51.555228},
      {50.358616, 50.358616, 62.999448, 62.999448, 50.358616, 50.358616},
      {49.071641, 49.071641, 62.999349, 62.999349, 49.071641, 49.071641},
      {47.669249, 47.669249, 62.999238, 62.999238, 47.669249, 47.669249},
      {46.111621, 46.111621, 62.999115, 62.999115, 46.111621, 46.111621},
      {44.327683, 44.327683, 62.998976, 62.998976, 44.327683, 44.327683},
      {42.461213, 42.461213, 62.998823, 62.998823, 42.461213, 42.461213},
      {40.554755, 40.554755, 62.998655, 62.998655, 40.554755, 40.554755},
      {38.605485, 38.605485, 62.998472, 62.998472, 38.605485, 38.605485},
      {36.610228, 36.610228, 62.998274, 62.998274, 36.610228, 36.610228},
      {34.565397, 34.565397, 62.998059, 62.998059, 34.565397, 34.565397},
      {32.466909, 32.466909, 62.997828, 62.997828, 32.466909, 32.466909}};

    EXPECT_EQ(waterContentSolution.size(), correctWaterContentSolution.size());

    for(auto i = 0u; i < waterContentSolution.size(); ++i)
    {
        for(auto j = 0u; j < waterContentSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctWaterContentSolution[i][j], waterContentSolution[i][j], 1e-6);
        }
    }

    std::vector<std::vector<double>> correctTemperatureSolution{
      {0.104018, 0.104018, 0.071546, 0.071546, 0.104018, 0.104018},
      {0.175125, 0.175125, 0.142795, 0.142795, 0.175125, 0.175125},
      {0.231102, 0.231102, 0.203541, 0.203541, 0.231102, 0.231102},
      {0.27688, 0.27688, 0.25399, 0.25399, 0.27688, 0.27688},
      {0.314786, 0.314786, 0.29581, 0.29581, 0.314786, 0.314786},
      {0.346416, 0.346416, 0.330618, 0.330618, 0.346416, 0.346416},
      {0.373016, 0.373016, 0.359776, 0.359776, 0.373016, 0.373016},
      {0.395589, 0.395589, 0.384403, 0.384403, 0.395589, 0.395589},
      {0.414954, 0.414954, 0.405408, 0.405408, 0.414954, 0.414954},
      {0.431777, 0.431777, 0.423536, 0.423536, 0.431777, 0.431777},
      {0.446607, 0.446607, 0.439395, 0.439395, 0.446607, 0.446607},
      {0.459896, 0.459896, 0.453488, 0.453488, 0.459896, 0.459896},
      {0.472024, 0.472024, 0.466231, 0.466231, 0.472024, 0.472024},
      {0.483313, 0.483313, 0.477979, 0.477979, 0.483313, 0.483313},
      {0.494048, 0.494048, 0.489038, 0.489038, 0.494048, 0.494048},
      {0.504496, 0.504496, 0.499687, 0.499687, 0.504496, 0.504496},
      {0.514929, 0.514929, 0.510202, 0.510202, 0.514929, 0.514929},
      {0.525667, 0.525667, 0.52089, 0.52089, 0.525667, 0.525667},
      {0.536711, 0.536711, 0.531837, 0.531837, 0.536711, 0.536711},
      {0.548039, 0.548039, 0.543061, 0.543061, 0.548039, 0.548039},
      {0.559656, 0.559656, 0.554572, 0.554572, 0.559656, 0.559656},
      {0.571573, 0.571573, 0.566381, 0.566381, 0.571573, 0.571573},
      {0.583806, 0.583806, 0.578503, 0.578503, 0.583806, 0.583806},
      {0.596373, 0.596373, 0.590954, 0.590954, 0.596373, 0.596373}};

    EXPECT_EQ(temperatureSolution.size(), correctTemperatureSolution.size());

    for(auto i = 0u; i < temperatureSolution.size(); ++i)
    {
        for(auto j = 0u; j < temperatureSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctTemperatureSolution[i][j], temperatureSolution[i][j], 1e-6);
        }
    }
}
