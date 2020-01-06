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

    // Material Properties (Cottaer Sandstone)
    auto & material =
      MaterialPool::Instance().createSolidMaterial("Test material");

    HygroThermFEM::MultiDomain domain;

    const auto matCheck = domain.checkMaterialsForTransientSimulation();
    EXPECT_EQ(matCheck.size(), 1);
}
