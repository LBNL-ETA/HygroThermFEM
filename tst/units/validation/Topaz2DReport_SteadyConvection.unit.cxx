#include <gtest/gtest.h>

#include "DumpCsv.hxx"
#include "HygroThermFEM2D.hxx"
#include "SlabCreator.hxx"
#include "TestMaterials.hxx"

/////////////////////////////////////////////////////////////////////////////////////
/// TOPAZ2D Validation Status Report (Drayton, UCRL-ID-106550, December 1990),
/// section 1.3 -- Plane Steady Linear, Convection Boundary Condition.
///
/// Unit slab (k = 1, L = 1) with the surface at x = 0 kept at 1.0 and convection
/// h (T - 0) with h = 1 at x = 1. The steady profile is linear,
/// T(x) = 1 - x / 2, and linear elements reproduce it nodally exactly.
/////////////////////////////////////////////////////////////////////////////////////

TEST(Topaz2DReport_SteadyConvection, LinearProfile)
{
    SCOPED_TRACE("Begin Test: steady slab, fixed temperature left, convection right.");

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

    constexpr auto tSurface = 1.0;
    multiDomain.thermal().createBC_FixedTemperature(1, 2, tSurface);

    constexpr auto tAir = 0.0;
    constexpr auto hc = 1.0;
    const HygroThermFEM::FixedBCHCCoefficients bcCoeff{tAir, hc};
    multiDomain.thermal().createBC_FixedHc(21, 22, bcCoeff);

    const auto solution = multiDomain.steadyState();

    TestHelper::CsvDump dump("topaz_steady_convection.csv", 11);
    dump.addRow(1, TestHelper::bottomRow(solution.temperature, 11, 2));

    // T(x) = Ts1 + h (Ta - Ts1) / (k + h L) * x = 1 - x / 2, exact for linear elements.
    for(std::size_t col = 0; col < 11; ++col)
    {
        const double expected = 1.0 - 0.05 * static_cast<double>(col);
        EXPECT_NEAR(expected, solution.temperature[col * 2], 1e-6) << "column " << col;
    }
}
