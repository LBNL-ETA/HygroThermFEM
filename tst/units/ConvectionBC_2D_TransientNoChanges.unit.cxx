#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

class ConvectionBC_2D_TransientNoChanges : public testing::Test
{
protected:
    void SetUp() override
    {}

    void TearDown() override
    {}
};

TEST_F(ConvectionBC_2D_TransientNoChanges, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Three elements with simple convection BC.");

    HygroThermFEM::MultiDomain multiDomain(true, false);

    std::vector<double> gridXCoordinates{0, 0.05, 0.1, 0.15};

    const double initialTemperature = 20;
    const double initialMoistureContent = 0;
    const double initialPressure = 101325;

    HygroThermFEM::State state(initialTemperature, initialMoistureContent, initialPressure, 0);
    size_t nodeIndex = 0;
    for(auto val : gridXCoordinates)
    {
        ++nodeIndex;
        multiDomain.nodes().createNode(nodeIndex, val, 0.00, state);
        ++nodeIndex;
        multiDomain.nodes().createNode(nodeIndex, val, 0.05, state);
    }

    // Material Properties (Cottaer Sandstone)
    constexpr double thermalConductivityDry{1.8};
    constexpr double density{2050.0};
    constexpr double porosity{0.22};
    constexpr double specificHeatCapacityDry{850.0};
    constexpr double diffusionResistanceFactor{15.0};
    const std::vector<FenestrationCommon::point> thermalConductivityMoistureDependent = {
      {0.0, 1.8}, {180, 1.8}};
    constexpr double thermalConductivityMeasuredAtTemperature{0};
    const std::vector<FenestrationCommon::point> thermalConductivityTemperatureDependent = {
      {0.0, 1.8}, {1, 1.8}};
    constexpr double thermalConductivityMeasuredAtHumidity{0};
    const std::vector<FenestrationCommon::point> liquidTransportationCurve = {{0, 0},
                                                                              {27, 1E-8},
                                                                              {45, 1.1E-8},
                                                                              {90, 2E-8},
                                                                              {126, 3.5E-8},
                                                                              {144, 5E-8},
                                                                              {162, 1E-7},
                                                                              {171, 2E-7},
                                                                              {180, 7E-7}};

    const std::vector<FenestrationCommon::point> moistureStorageFunction = {{0, 0},
                                                                            {0.5, 5.3},
                                                                            {0.65, 8.4},
                                                                            {0.8, 12},
                                                                            {0.93, 17},
                                                                            {0.95, 25},
                                                                            {0.99, 63},
                                                                            {0.995, 83},
                                                                            {0.999, 120},
                                                                            {1, 180}};

    auto & material =
      multiDomain.materials().createSolidMaterial("Cottaer Sandstone",
                                                  thermalConductivityDry,
                                                  density,
                                                  porosity,
                                                  specificHeatCapacityDry,
                                                  diffusionResistanceFactor,
                                                  thermalConductivityMoistureDependent,
                                                  thermalConductivityMeasuredAtTemperature,
                                                  thermalConductivityTemperatureDependent,
                                                  thermalConductivityMeasuredAtHumidity,
                                                  liquidTransportationCurve,
                                                  moistureStorageFunction);

    /// Create elements
    for(size_t idx = 1u; idx <= (multiDomain.nodes().maxIndex() - 2) / 2; ++idx)
    {
        const auto index1 = 2u * idx + 1u;
        const auto index2 = 2u * idx + 2u;
        const auto index3 = 2u * idx;
        const auto index4 = 2u * idx - 1u;
        multiDomain.createElement(index1, index2, index3, index4, material.name());
    }

    // Create Boundary Conditions
    constexpr auto tSurface = 20.0;
    constexpr auto hc = 1.0;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff{tSurface, hc};

    multiDomain.thermal().createBC_FixedHc(1, 2, bcCoeff);

    constexpr auto dTime = 36000;
    constexpr auto nSteps = 4;

    auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
    std::vector<std::vector<double>> solution;

    for(unsigned i = 0; i < nSteps; ++i)
    {
        temperatures = multiDomain.thermal().transient(temperatures, dTime).solution;
        solution.push_back(temperatures);
    }

    std::vector<std::vector<double>> correctSolution{
      {20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0},
      {20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0},
      {20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0},
      {20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0}};

    EXPECT_EQ(solution.size(), correctSolution.size());

    for(auto i = 0u; i < correctSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctSolution[i][j], solution[i][j], 1e-6);
        }
    }
}
