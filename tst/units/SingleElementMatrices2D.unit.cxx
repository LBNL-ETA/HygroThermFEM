#include <gtest/gtest.h>
#include <stdexcept>

#include "HygroThermFEM2D.hxx"

class TestSingleElementMatrices2D : public testing::Test
{
protected:
    void SetUp() override
    {}

    void TearDown() override
    {}
};

TEST_F(TestSingleElementMatrices2D, TestConductionMatrix)
{
    SCOPED_TRACE("Begin Test: Single element isothropic conduction matrix and "
                 "RhoCp matrix.");

    // Create MultiDomain before nodes
    HygroThermFEM::MultiDomain multiDomain(true, false);

    // Enter nodes. Arguments are: node number, x-coordinate, y-coordinate

    multiDomain.nodes().createNode(1, 5, 5);
    multiDomain.nodes().createNode(2, 5, 0);
    multiDomain.nodes().createNode(3, 15, 0);
    multiDomain.nodes().createNode(4, 15, 5);

    // Material Properties
    constexpr double thermalConductivityDry{1.0};
    constexpr double density{1.0};
    constexpr double porosity{0.0};
    constexpr double specificHeatCapacityDry{1.0};
    constexpr double diffusionResistanceFactor{15.0};
    const std::vector<FenestrationCommon::point> thermalConductivityMoistureDependent = {
      {0.0, 1.0}, {180, 1.0}};
    constexpr double thermalConductivityMeasuredAtTemperature{0};
    const std::vector<FenestrationCommon::point> thermalConductivityTemperatureDependent = {
      {0.0, 1.0}, {1, 1.0}};
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
      multiDomain.materials().createSolidMaterial("Test Material",
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

    const HygroThermFEM::ElementThermalLinear2D aElem{multiDomain.nodes(), multiDomain.materials(), 1, 2, 3, 4, material.name()};

    auto condMat = aElem.DDuMatrices();

    std::vector<std::vector<double>> correctCondMat = {
      {0.833333333, -0.583333333, -0.416666667, 0.166666667},
      {-0.583333333, 0.833333333, 0.166666667, -0.416666667},
      {-0.416666667, 0.166666667, 0.833333333, -0.583333333},
      {0.166666667, -0.416666667, -0.583333333, 0.833333333}};

    for(auto i = 0; i < 4; ++i)
    {
        for(auto j = 0; j < 4; ++j)
        {
            EXPECT_NEAR(correctCondMat[i][j], condMat(i, j), 1e-6);
        }
    }

    auto rhoCpMat = aElem.capacitanceMatrices();

    std::vector<std::vector<double>> correctRhoCpMat = {
      {5.555555556, 2.777777778, 1.388888889, 2.777777778},
      {2.777777778, 5.555555556, 2.777777778, 1.388888889},
      {1.388888889, 2.777777778, 5.555555556, 2.777777778},
      {2.777777778, 1.388888889, 2.777777778, 5.555555556}};

    for(auto i = 0; i < 4; ++i)
    {
        for(auto j = 0; j < 4; ++j)
        {
            EXPECT_NEAR(correctRhoCpMat[i][j], rhoCpMat(i, j), 1e-6);
        }
    }
}
