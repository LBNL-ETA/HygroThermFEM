#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::NodePool;
using HygroThermFEM::MaterialPool;
using HygroThermFEM::State;

class TestMaterialChecker : public testing::Test
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

TEST_F(TestMaterialChecker, TestExample_1)
{
    SCOPED_TRACE("Test for material missing properties.");

    // Note that material is incomplete
    MaterialPool::Instance().createSolidMaterial("Test material");

    HygroThermFEM::MultiDomain domain;

    const auto matCheck = domain.checkMaterialsForTransientSimulation();
    EXPECT_EQ(matCheck.size(), 1);

    // Series of test means following. Material property is missing for given simulation properties
    // and therefore it is set to true. If material property is not needed for given simulation
    // properties then it is set to false.
    EXPECT_EQ(matCheck[0].Density, true);
    EXPECT_EQ(matCheck[0].Emissivity, true);
    EXPECT_EQ(matCheck[0].Porosity, true);
    EXPECT_EQ(matCheck[0].SpecificHeatCapacityDry, true);
    EXPECT_EQ(matCheck[0].ThermalConductivityDry, true);
    EXPECT_EQ(matCheck[0].WaterVaporDiffusionResistanceFactor, true);
    EXPECT_EQ(matCheck[0].MoistureStorageFunction, true);
    EXPECT_EQ(matCheck[0].LiquidTransportationSuction, true);
    EXPECT_EQ(matCheck[0].LiquidTransportationRedistribution, false);
    EXPECT_EQ(matCheck[0].ThermalConductivityMoistureAndTemperatureDependent, false);
}
