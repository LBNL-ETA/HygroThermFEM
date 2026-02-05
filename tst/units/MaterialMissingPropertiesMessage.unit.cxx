#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::MaterialMissingProperties;

TEST(TestMaterialMissingPropertiesMessage, NoMissingProperties)
{
    SCOPED_TRACE("Begin Test: No missing properties returns empty message.");

    MaterialMissingProperties props;
    props.materialName = "TestMaterial";
    // All flags default to false

    const auto message = props.missingPropertiesMessage();

    EXPECT_TRUE(message.empty());
}

TEST(TestMaterialMissingPropertiesMessage, MissingDensity)
{
    SCOPED_TRACE("Begin Test: Missing Density property.");

    MaterialMissingProperties props;
    props.materialName = "TestMaterial";
    props.Density = true;

    const auto message = props.missingPropertiesMessage();

    ASSERT_EQ(message.size(), 2u);
    EXPECT_EQ(message[0], "Material TestMaterial is missing following properties:");
    EXPECT_EQ(message[1], "- Density");
}

TEST(TestMaterialMissingPropertiesMessage, MissingEmissivity)
{
    SCOPED_TRACE("Begin Test: Missing Emissivity property.");

    MaterialMissingProperties props;
    props.materialName = "TestMaterial";
    props.Emissivity = true;

    const auto message = props.missingPropertiesMessage();

    ASSERT_EQ(message.size(), 2u);
    EXPECT_EQ(message[1], "- Emissivity");
}

TEST(TestMaterialMissingPropertiesMessage, MissingPorosity)
{
    SCOPED_TRACE("Begin Test: Missing Porosity property.");

    MaterialMissingProperties props;
    props.materialName = "TestMaterial";
    props.Porosity = true;

    const auto message = props.missingPropertiesMessage();

    ASSERT_EQ(message.size(), 2u);
    EXPECT_EQ(message[1], "- Porosity");
}

TEST(TestMaterialMissingPropertiesMessage, MissingSpecificHeatCapacityDry)
{
    SCOPED_TRACE("Begin Test: Missing SpecificHeatCapacityDry property.");

    MaterialMissingProperties props;
    props.materialName = "TestMaterial";
    props.SpecificHeatCapacityDry = true;

    const auto message = props.missingPropertiesMessage();

    ASSERT_EQ(message.size(), 2u);
    EXPECT_EQ(message[1], "- Specific Heat Capacity Dry");
}

TEST(TestMaterialMissingPropertiesMessage, MissingThermalConductivityDry)
{
    SCOPED_TRACE("Begin Test: Missing ThermalConductivityDry property.");

    MaterialMissingProperties props;
    props.materialName = "TestMaterial";
    props.ThermalConductivityDry = true;

    const auto message = props.missingPropertiesMessage();

    ASSERT_EQ(message.size(), 2u);
    EXPECT_EQ(message[1], "- Thermal Conductivity Dry");
}

TEST(TestMaterialMissingPropertiesMessage, MissingWaterVaporDiffusionResistanceFactor)
{
    SCOPED_TRACE("Begin Test: Missing WaterVaporDiffusionResistanceFactor property.");

    MaterialMissingProperties props;
    props.materialName = "TestMaterial";
    props.WaterVaporDiffusionResistanceFactor = true;

    const auto message = props.missingPropertiesMessage();

    ASSERT_EQ(message.size(), 2u);
    EXPECT_EQ(message[1], "- Water Vapor Diffusion Resistance Factor");
}

TEST(TestMaterialMissingPropertiesMessage, MissingMoistureStorageFunction)
{
    SCOPED_TRACE("Begin Test: Missing MoistureStorageFunction property.");

    MaterialMissingProperties props;
    props.materialName = "TestMaterial";
    props.MoistureStorageFunction = true;

    const auto message = props.missingPropertiesMessage();

    ASSERT_EQ(message.size(), 2u);
    EXPECT_EQ(message[1], "- Moisture Storage Function");
}

TEST(TestMaterialMissingPropertiesMessage, MissingLiquidTransportationSuction)
{
    SCOPED_TRACE("Begin Test: Missing LiquidTransportationSuction property.");

    MaterialMissingProperties props;
    props.materialName = "TestMaterial";
    props.LiquidTransportationSuction = true;

    const auto message = props.missingPropertiesMessage();

    ASSERT_EQ(message.size(), 2u);
    EXPECT_EQ(message[1], "- Liquid Transportation Suction Curve");
}

TEST(TestMaterialMissingPropertiesMessage, MissingLiquidTransportationRedistribution)
{
    SCOPED_TRACE("Begin Test: Missing LiquidTransportationRedistribution property.");

    MaterialMissingProperties props;
    props.materialName = "TestMaterial";
    props.LiquidTransportationRedistribution = true;

    const auto message = props.missingPropertiesMessage();

    ASSERT_EQ(message.size(), 2u);
    EXPECT_EQ(message[1], "- Liquid Transportation Redistribution Curve");
}

TEST(TestMaterialMissingPropertiesMessage, MissingThermalConductivityMoistureAndTemperatureDependent)
{
    SCOPED_TRACE("Begin Test: Missing ThermalConductivityMoistureAndTemperatureDependent property.");

    MaterialMissingProperties props;
    props.materialName = "TestMaterial";
    props.ThermalConductivityMoistureAndTemperatureDependent = true;

    const auto message = props.missingPropertiesMessage();

    ASSERT_EQ(message.size(), 2u);
    EXPECT_EQ(message[1], "- Thermal Conductivity Moisture and Temperature Dependent");
}

TEST(TestMaterialMissingPropertiesMessage, MultiplePropertiesMissing)
{
    SCOPED_TRACE("Begin Test: Multiple properties missing.");

    MaterialMissingProperties props;
    props.materialName = "ConcreteMaterial";
    props.Density = true;
    props.Porosity = true;
    props.MoistureStorageFunction = true;

    const auto message = props.missingPropertiesMessage();

    ASSERT_EQ(message.size(), 4u);
    EXPECT_EQ(message[0], "Material ConcreteMaterial is missing following properties:");
    EXPECT_EQ(message[1], "- Density");
    EXPECT_EQ(message[2], "- Porosity");
    EXPECT_EQ(message[3], "- Moisture Storage Function");
}

TEST(TestMaterialMissingPropertiesMessage, AllPropertiesMissing)
{
    SCOPED_TRACE("Begin Test: All properties missing.");

    MaterialMissingProperties props;
    props.materialName = "EmptyMaterial";
    props.Density = true;
    props.Emissivity = true;
    props.Porosity = true;
    props.SpecificHeatCapacityDry = true;
    props.ThermalConductivityDry = true;
    props.WaterVaporDiffusionResistanceFactor = true;
    props.MoistureStorageFunction = true;
    props.LiquidTransportationSuction = true;
    props.LiquidTransportationRedistribution = true;
    props.ThermalConductivityMoistureAndTemperatureDependent = true;

    const auto message = props.missingPropertiesMessage();

    // Header + 10 properties
    ASSERT_EQ(message.size(), 11u);
    EXPECT_EQ(message[0], "Material EmptyMaterial is missing following properties:");
}
