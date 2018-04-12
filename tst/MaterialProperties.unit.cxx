#include <memory>
#include <gtest/gtest.h>

#include "MoisThermFEM2D.hxx"

using MoisThermFEM::Property;
using MoisThermFEM::State;
using MoisThermFEM::MaterialProperties;

class MaterialPropertiesUnit : public testing::Test {

protected:
	void
	SetUp() override {
	}

	void
	TearDown() override {
		MoisThermFEM::MaterialPool::Instance().clear();
	}

};

TEST_F( MaterialPropertiesUnit, TestMaterialProperties ) {
	SCOPED_TRACE( "Begin Test: Test water and moisture fill." );

	auto & materialPool = MoisThermFEM::MaterialPool::Instance();

	auto & material = materialPool.Instance().createMaterial(
			"Cottaer Sandstone",
			2050, /// density
			0.22, /// porosity
			850,  /// specific heat capacity (dry)
			1.8,  /// thermal conductivity (dry)
			15,   /// diffusion resistance factor
			{ { 0,   0 },  /// liquid transportation coefficient
				{ 27,  1E-8 },
				{ 45,  1.1E-8 },
				{ 90,  2E-8 },
				{ 126, 3.5E-8 },
				{ 144, 5E-8 },
				{ 162, 1E-7 },
				{ 171, 2E-7 },
				{ 180, 7E-7 } },
			{ { 0,     0 },   /// sorption curve
				{ 0.5,   5.3 },
				{ 0.65,  8.4 },
				{ 0.8,   12 },
				{ 0.93,  17 },
				{ 0.95,  25 },
				{ 0.99,  63 },
				{ 0.995, 83 },
				{ 0.999, 120 },
				{ 1,     180 } }
	);

	auto airFill = MaterialProperties::getAirFill( material );
	auto waterFill = MaterialProperties::getWaterFill( material );

	State state1( 273.15, 0.0, 101325 );

	auto aFill = airFill->value( state1 );
	auto wFill = waterFill->value( state1 );

	EXPECT_NEAR( 0.22, aFill, 1e-6 );
	EXPECT_NEAR( 0.00, wFill, 1e-6 );

	State state2( 273.15, 1.0, 101325 );

	aFill = airFill->value( state2 );
	wFill = waterFill->value( state2 );

	EXPECT_NEAR( 0.00, aFill, 1e-6 );
	EXPECT_NEAR( 0.22, wFill, 1e-6 );

	State state3( 273.15, 0.5, 101325 );

	aFill = airFill->value( state3 );
	wFill = waterFill->value( state3 );

	EXPECT_NEAR( 0.213522, aFill, 1e-6 );
	EXPECT_NEAR( 0.006478, wFill, 1e-6 );

}
