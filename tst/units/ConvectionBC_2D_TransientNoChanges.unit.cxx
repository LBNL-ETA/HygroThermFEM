#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::NodePool;
using HygroThermFEM::MaterialPool;

class ConvectionBC_2D_TransientNoChanges : public testing::Test
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

TEST_F(ConvectionBC_2D_TransientNoChanges, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Three elements with simple convection BC.");

    std::vector<double> gridXCoordinates{0, 0.05, 0.1, 0.15};

    const double initialTemperature = 20;
    const double initialMoistureContent = 0;
    const double initialPressure = 101325;

    HygroThermFEM::State state(initialTemperature, initialMoistureContent, initialPressure, 0);
    size_t nodeIndex = 0;
    for(auto val : gridXCoordinates)
    {
        ++nodeIndex;
        NodePool::Instance().createNode(nodeIndex, val, 0.00, state);
        ++nodeIndex;
        NodePool::Instance().createNode(nodeIndex, val, 0.05, state);
    }

    // Material Properties (Cottaer Sandstone)
    const double thermalConductivityDry{1.8};
    const double density{2050.0};
    const double porosity{0.22};
    const double specificHeatCapacityDry{850.0};
    const double diffusionResistanceFactor{15.0};
    const std::vector<FenestrationCommon::point> thermalConductivityMoistureDependent = {
      {0.0, 1.8}, {180, 1.8}};
    const double thermalConductivityMeasuredAtTemperature{0};
    const std::vector<FenestrationCommon::point> thermalConductivityTemperatureDependent = {
      {0.0, 1.8}, {1, 1.8}};
    const double thermalConductivityMeasuredAtHumidity{0};
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

    HygroThermFEM::SingleDomain domain{HygroThermFEM::DomainType::Thermal};

    /// Create elements
    for(size_t i = 1u; i <= (NodePool::Instance().maxIndex() - 2) / 2; ++i)
    {
        const auto index1 = 2u * i + 1u;
        const auto index2 = 2u * i + 2u;
        const auto index3 = 2u * i;
        const auto index4 = 2u * i - 1u;
        createElement(domain, index1, index2, index3, index4, material.name());
    }

    // Create Boundary Conditions
    const auto tSurface = 20.0;
    const auto hc = 1.0;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff{tSurface, hc};

    HygroThermFEM::Thermal::createBC_FixedHc(domain, 1, 2, bcCoeff);

    const auto dTime = 36000;
    const auto nSteps = 4;

    auto temperatures = properties(HygroThermFEM::Variable::temperature);
    std::vector<std::vector<double>> solution;

    for(unsigned i = 0; i < nSteps; ++i)
    {
        temperatures = HygroThermFEM::Substitution::transient(domain, temperatures, dTime).solution;
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
