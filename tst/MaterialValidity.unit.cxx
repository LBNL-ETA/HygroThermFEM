#include <gtest/gtest.h>
#include <memory>
#include <vector>

#include "HygroThermFEM2D.hxx"

class MaterialValidityTest : public testing::Test
{
protected:
    void SetUp() override
    {
        HygroThermFEM::MaterialPool::Instance().createSolidMaterial(materialName());
    }

    void TearDown() override
    {
        HygroThermFEM::NodePool::Instance().clear();
        HygroThermFEM::MaterialPool::Instance().clear();
    }

    std::string materialName() const
    {
        return m_MaterialName;
    }

private:
    std::string m_MaterialName{"Test"};
};

TEST_F(MaterialValidityTest, TestConductivity)
{
    SCOPED_TRACE("Test for material conductivity.");
    const auto & mat = HygroThermFEM::MaterialPool::Instance().material(materialName());
    try
    {
        mat.thermalConductivityDry();
    }
    catch(const std::exception & e)
    {
        // and this tests that it has the correct message
        EXPECT_EQ("Material " + materialName() + " do not have thermal conductivity.", e.what());
    }
}

TEST_F(MaterialValidityTest, TestDensity)
{
    SCOPED_TRACE("Test for material conductivity.");
    const auto & mat = HygroThermFEM::MaterialPool::Instance().material(materialName());
    try
    {
        mat.density();
    }
    catch(const std::exception & e)
    {
        // and this tests that it has the correct message
        EXPECT_EQ("Material " + materialName() + " do not have density.", e.what());
    }
}

TEST_F(MaterialValidityTest, TestHeatCapacity)
{
    SCOPED_TRACE("Test for material specific heat capacity.");
    const auto & mat = HygroThermFEM::MaterialPool::Instance().material(materialName());
    try
    {
        mat.heatCapacity();
    }
    catch(const std::exception & e)
    {
        // and this tests that it has the correct message
        EXPECT_EQ("Material " + materialName() + " do not have specific heat capacity.", e.what());
    }
}

TEST_F(MaterialValidityTest, TestPorosity)
{
    SCOPED_TRACE("Test for material porosity.");
    const auto & mat = HygroThermFEM::MaterialPool::Instance().material(materialName());
    try
    {
        mat.porosity();
    }
    catch(const std::exception & e)
    {
        // and this tests that it has the correct message
        EXPECT_EQ("Material " + materialName() + " do not have porosity.", e.what());
    }
}

TEST_F(MaterialValidityTest, TestEmissivity)
{
    SCOPED_TRACE("Test for material emissivity.");
    const auto & mat = HygroThermFEM::MaterialPool::Instance().material(materialName());
    try
    {
        mat.emissivity();
    }
    catch(const std::exception & e)
    {
        // and this tests that it has the correct message
        EXPECT_EQ("Material " + materialName() + " do not have emissivity.", e.what());
    }
}

TEST_F(MaterialValidityTest, TestDiffusionResistanceFactor)
{
    SCOPED_TRACE("Test for material diffusion resistance factor.");
    const auto & mat = HygroThermFEM::MaterialPool::Instance().material(materialName());
    try
    {
        mat.diffusionResistanceFactor();
    }
    catch(const std::exception & e)
    {
        // and this tests that it has the correct message
        EXPECT_EQ("Material " + materialName() + " do not have diffusion resistance factor.", e.what());
    }
}