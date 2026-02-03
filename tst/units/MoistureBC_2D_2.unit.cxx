#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::NodePool;
using HygroThermFEM::MaterialPool;

class MoistureBC_2D_2 : public testing::Test
{
protected:
    void SetUp() override
    {}

    void TearDown() override
    {
        NodePool::Instance().clear();
        MaterialPool::Instance().clear();
    }
};

TEST_F(MoistureBC_2D_2, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Simple two elements example with moisture transfer.");

    std::vector<double> gridXCoordinates{0.15, 0.05, 0.00};

    constexpr auto domainTemperature = 20.0;
    constexpr auto domainHumidity = 0.0;
    constexpr auto domainPressure = 101325;
    constexpr auto liquidPercentage = 1.0;

    const HygroThermFEM::State state(
      domainTemperature, domainHumidity, domainPressure, liquidPercentage);
    size_t nodeIndex = 0;
    for(auto val : gridXCoordinates)
    {
        ++nodeIndex;
        NodePool::Instance().createNode(nodeIndex, val, 0.05, state);
        ++nodeIndex;
        NodePool::Instance().createNode(nodeIndex, val, 0.00, state);
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
      MaterialPool::Instance().createSolidMaterial("Cottaer Sandstone",
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

    HygroThermFEM::MultiDomain multiDomain(false, true);
    auto & domain = multiDomain.moisture();

    /// Create elements
    for(size_t i = 1; i <= (NodePool::Instance().maxIndex() - 2) / 2; ++i)
    {
        const auto node1 = 2u * i + 1u;
        const auto node2 = 2u * i + 2u;
        const auto node3 = 2u * i;
        const auto node4 = 2u * i - 1u;
        multiDomain.createElement(node2, node3, node4, node1, material.name());
    }

    // Create Boundary Conditions
    constexpr auto ambientTemperature = 20.0;
    constexpr auto ambientHumidity = 0.2;
    const auto surfaceTilt{90.0};

    const HygroThermFEM::TARPCoefficients bcCoeff{ambientTemperature, ambientHumidity};

    domain.createBC_TARPHc(5, 6, bcCoeff, surfaceTilt);

    constexpr auto dTime = 3600;
    constexpr auto nSteps = 4;

    auto humidities = NodePool::Instance().properties(HygroThermFEM::Variable::humidity);
    std::vector<std::vector<double>> solution;

    for(size_t i = 0u; i < nSteps; ++i)
    {
        humidities = domain.transient(humidities, dTime).solution;
        auto waterContent = NodePool::Instance().properties(HygroThermFEM::Variable::water);
        solution.push_back(waterContent);
    }

    std::vector<std::vector<double>> correctSolution{
      {3.748586e-06, 3.748586e-06, 0.001924969331, 0.001924969331, 0.7428227662, 0.7428227662},
      {1.364027149e-05, 1.364027149e-05, 0.005083314705, 0.005083314705, 1.221651619, 1.221651619},
      {3.11968812e-05, 3.11968812e-05, 0.009029288908, 0.009029288908, 1.530317721, 1.530317721},
      {5.737783803e-05, 5.737783803e-05, 0.01347560965, 0.01347560965, 1.729301858, 1.729301858}};

    EXPECT_EQ(solution.size(), correctSolution.size());

    for(auto i = 0u; i < correctSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctSolution[i][j], solution[i][j], 1e-9);
        }
    }
}
