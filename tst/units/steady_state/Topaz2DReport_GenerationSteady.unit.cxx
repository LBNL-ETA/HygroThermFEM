#include <gtest/gtest.h>

#include "DumpCsv.hxx"
#include "HygroThermFEM2D.hxx"
#include "SlabCreator.hxx"
#include "TestMaterials.hxx"

/////////////////////////////////////////////////////////////////////////////////////
/// TOPAZ2D Validation Status Report (Drayton, UCRL-ID-106550, December 1990),
/// section 1.4 -- Plane Steady Linear, Internal Energy Generation.
///
/// Unit slab (k = 1, L = 1) with uniform volumetric generation q = 1, adiabatic at
/// x = 0, surface at x = 1 held at 0. The steady profile is the parabola
///
///     T(x) = (q L^2 / 2k) (1 - (x/L)^2)
///
/// (Incropera & DeWitt p. 86). With the consistent element load and a uniform mesh
/// the nodal solution reproduces the parabola exactly: the discrete second
/// difference of a quadratic is exact.
/////////////////////////////////////////////////////////////////////////////////////

TEST(Topaz2DReport_GenerationSteady, ParabolicProfile)
{
    SCOPED_TRACE("Begin Test: steady slab with uniform internal generation.");

    HygroThermFEM::MultiDomain multiDomain({.performThermal = true, .performMoisture = false});

    constexpr HygroThermFEM::State state({
        .temperature = 0.0,
        .humidity = 0.0,
        .pressure = 101325.0,
        .liquidPercent = 0
    });

    auto params = TestHelper::TestMaterial();
    params.density = 1.0;
    params.heatCapacity = 1.0;
    const auto & material = multiDomain.materials().createSolidMaterial(params);

    TestHelper::SlabBuilder(multiDomain)
        .gridXCoordinates({0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0})
        .height(0.05)
        .material(material.name())
        .state(state)
        .build();

    constexpr auto tSurface = 0.0;
    multiDomain.thermal().createBC_FixedTemperature(21, 22, tSurface);

    constexpr auto generation = 1.0;
    multiDomain.thermal().setVolumetricSource(generation);

    const auto solution = multiDomain.steadyState();

    TestHelper::CsvDump dump("topaz_generation_steady.csv", 11);
    dump.addRow(1, TestHelper::bottomRow(solution.temperature, 11, 2));

    for(std::size_t col = 0; col < 11; ++col)
    {
        const double position = 0.1 * static_cast<double>(col);
        const double expected = 0.5 * (1.0 - position * position);
        EXPECT_NEAR(expected, solution.temperature[col * 2], 1e-6) << "column " << col;
    }
}
