#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "SlabCreator.hxx"
#include "TestMaterials.hxx"

/////////////////////////////////////////////////////////////////////////////////////
/// TOPAZ2D Validation Status Report (Drayton, UCRL-ID-106550, December 1990),
/// section 2.3 -- Plane Steady Nonlinear, Radiation Boundary Condition.
///
/// Unit slab (k = 1, L = 1) with the surface at x = 0 kept at 0 C and black-body
/// radiative exchange (eps = 0.5) with a 100 C environment at x = 1. The steady
/// profile is linear, T(x) = C0 x, with C0 the root of
///
///     eps sigma ((C0 L + 273.15)^4 - 373.15^4) + k C0 = 0
///
/// The report states the same problem with temperatures already absolute; here the
/// parameters are engine-native (Celsius, sigma = 5.6697e-8), and the expected
/// slope C0 = 84.707885049305 comes from the same root equation in those units
/// (hygrothermfem_python, analytic.radiation_steady_slope). The steady coupling
/// loop re-linearizes the radiation coefficient about each pass's surface
/// temperature, which is the Picard iteration for the T^4 balance.
/////////////////////////////////////////////////////////////////////////////////////

TEST(Topaz2DReport_SteadyRadiation, LinearProfileFromRadiativeBalance)
{
    SCOPED_TRACE("Begin Test: steady slab, fixed temperature left, black-body right.");

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
    multiDomain.thermal().createBC_FixedTemperature(1, 2, tSurface);

    constexpr auto emissivity = 0.5;
    constexpr auto tRadiation = 100.0;
    multiDomain.thermal().createBC_BlackBodyRadiation(21, 22, emissivity, tRadiation);

    const auto solution = multiDomain.steadyState();


    constexpr double slope = 84.707885049305;
    for(std::size_t col = 0; col < 11; ++col)
    {
        const double expected = slope * 0.1 * static_cast<double>(col);
        EXPECT_NEAR(expected, solution.temperature[col * 2], 1e-3) << "column " << col;
    }
}
