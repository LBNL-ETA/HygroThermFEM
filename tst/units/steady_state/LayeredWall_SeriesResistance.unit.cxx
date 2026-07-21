#include <gtest/gtest.h>

#include "BeamBuilder.hxx"
#include "HygroThermFEM2D.hxx"
#include "TestMaterials.hxx"

/////////////////////////////////////////////////////////////////////////////////////
/// Two-layer wall between fixed surface temperatures: the steady profile is
/// piecewise linear with the interface temperature set by the series resistance
///
///     T_int = T_left + (T_right - T_left) * R_1 / (R_1 + R_2),  R_i = d_i / k_i
///
/// Cottaer Sandstone (k = 1.8, 0.06 m) against Stucco (k = 0.85, 0.04 m).
/// Nodally exact for linear elements, so the tolerance is a solver tolerance,
/// not a discretization band. The bottom node row is the layered 1D proxy the
/// validation book compares (hygrothermfem_python tests/test_multimaterial.py
/// asserts the same closed form for the reference solver).
/////////////////////////////////////////////////////////////////////////////////////

TEST(LayeredWall_SeriesResistance, TwoLayerSteadyConduction)
{
    SCOPED_TRACE("Begin Test: layered wall, fixed temperatures both faces.");

    HygroThermFEM::MultiDomain multiDomain({.performThermal = true, .performMoisture = false});

    constexpr HygroThermFEM::State state({
        .temperature = 0.0,
        .humidity = 0.0,
        .pressure = 101325.0,
        .liquidPercent = 0
    });

    const auto & cottaer =
      multiDomain.materials().createSolidMaterial(TestHelper::CottaerSandstone());
    const auto & stucco = multiDomain.materials().createSolidMaterial(TestHelper::Stucco());

    TestHelper::BeamBuilder builder(multiDomain);
    builder.xStart(0.0)
      .height(0.01)
      .numElementsY(1)
      .state(state)
      .addSegment({.material = cottaer.name(), .numElementsX = 6, .width = 0.06})
      .addSegment({.material = stucco.name(), .numElementsX = 4, .width = 0.04})
      .build();

    constexpr auto tLeft = 20.0;
    constexpr auto tRight = 0.0;
    for(const auto [index1, index2] : builder.leftEdge())
    {
        multiDomain.thermal().createBC_FixedTemperature(index1, index2, tLeft);
    }
    for(const auto [index1, index2] : builder.rightEdge())
    {
        multiDomain.thermal().createBC_FixedTemperature(index1, index2, tRight);
    }

    const auto solution = multiDomain.steadyState();


    // Series resistance: R1 = 0.06 / 1.8, R2 = 0.04 / 0.85.
    constexpr auto resistance1 = 0.06 / 1.8;
    constexpr auto resistance2 = 0.04 / 0.85;
    constexpr auto tInterface = tLeft + (tRight - tLeft) * resistance1
                                          / (resistance1 + resistance2);

    for(std::size_t col = 0; col < 11; ++col)
    {
        const double position = 0.01 * static_cast<double>(col);
        const double expected =
          position <= 0.06
            ? tLeft + (tInterface - tLeft) * position / 0.06
            : tInterface + (tRight - tInterface) * (position - 0.06) / 0.04;
        EXPECT_NEAR(expected, solution.temperature[col * 2], 1e-6) << "column " << col;
    }
}
