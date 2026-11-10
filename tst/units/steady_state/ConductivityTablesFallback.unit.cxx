#include <cmath>

#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

#include "BeamBuilder.hxx"
#include "TestMaterials.hxx"

namespace
{
    //! Stucco with its moisture and temperature conductivity tables removed: a material
    //! that carries only a dry conductivity, as a metal or a glazing layer does.
    HygroThermFEM::SolidMaterialParams stuccoWithoutTables()
    {
        auto params{TestHelper::Stucco()};
        params.thermalConductivityMoistureDependent = {};
        params.thermalConductivityTemperatureDependent = {};
        return params;
    }

    //! A beam of that material between a cold and a warm film, heat transfer alone.
    void buildBeam(HygroThermFEM::MultiDomain & multiDomain)
    {
        const auto & material = multiDomain.materials().createSolidMaterial(stuccoWithoutTables());
        TestHelper::BeamBuilder builder(multiDomain);
        builder.xStart(0.0)
          .height(0.05)
          .numElementsY(1)
          .state({.temperature = 10.0, .humidity = 0.5, .pressure = 101325.0, .liquidPercent = 1.0})
          .addSegment({.material = material.name(), .numElementsX = 10, .width = 0.2})
          .build();
        builder.applyBC_FixedHc(TestHelper::BeamBuilder::Edge::Left,
                                HygroThermFEM::FixedBCHCCoefficients{-18.0, 20.0, 0.0});
        builder.applyBC_FixedHc(TestHelper::BeamBuilder::Edge::Right,
                                HygroThermFEM::FixedBCHCCoefficients{21.0, 10.0, 0.0});
    }

    HygroThermFEM::Solution solve(const bool tablesAskedFor)
    {
        HygroThermFEM::MultiDomain multiDomain({.performThermal = true, .performMoisture = false});
        multiDomain.setPhysicsOptions({.thermalConductivityMoistureAndTemperatureDependent = tablesAskedFor});
        buildBeam(multiDomain);
        return multiDomain.steadyState();
    }
}   // namespace

// The conductivity tables are used where a material carries them. A material without
// them conducts at its dry value: it used to get no conduction term at all, and the
// solve failed to factorize the system matrix.
TEST(ConductivityTablesFallback, AMaterialWithoutTablesConductsAtItsDryValue)
{
    const auto withOption = solve(true);
    const auto withoutOption = solve(false);

    ASSERT_TRUE(withOption.converged);
    ASSERT_EQ(withOption.temperature.size(), withoutOption.temperature.size());
    for(std::size_t node = 0u; node < withOption.temperature.size(); ++node)
    {
        ASSERT_TRUE(std::isfinite(withOption.temperature[node]));
        EXPECT_NEAR(withOption.temperature[node], withoutOption.temperature[node], 1e-9) << "node " << node;
    }
}

// Such a material is complete for a heat transfer run with the option on: the tables are
// not required of it.
TEST(ConductivityTablesFallback, TheCheckerDoesNotAskForTables)
{
    HygroThermFEM::MultiDomain multiDomain({.performThermal = true, .performMoisture = false});
    multiDomain.setPhysicsOptions({.thermalConductivityMoistureAndTemperatureDependent = true});
    buildBeam(multiDomain);

    const auto check = multiDomain.checkForMaterialsValidity(HygroThermFEM::SimulationType::SteadyState);

    EXPECT_TRUE(check.isMaterialLibraryCorrect());
}
