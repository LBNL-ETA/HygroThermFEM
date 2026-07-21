#include <array>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "SlabCreator.hxx"
#include "TestMaterials.hxx"

/////////////////////////////////////////////////////////////////////////////////////
/// TOPAZ2D Validation Status Report (Drayton, UCRL-ID-106550, December 1990),
/// section 3.7 -- Plane Transient Linear, Temporal Convection Boundary Condition.
///
/// Unit slab (alpha = 1, k = 1) at zero initial temperature, adiabatic at x = 0,
/// convection h = 1 at x = 1 into an ambient that ramps as Ta = C0 t with C0 = 1.
/// Analytical solution: Carslaw & Jaeger p. 127. Expected values are the exact
/// series at the report's checkpoints (hygrothermfem_python,
/// analytic.slab_convection_ramp); the tolerance is the measured backward-Euler
/// discretization band at dt = 0.01 on ten elements (max 1.9e-3).
///
/// Per-step BC vectors advance only through the explicit timestepIndex argument
/// of MultiDomain::transient (element i = value during step i + 1), so the steps
/// are driven manually here rather than through transientMultiStep.
/////////////////////////////////////////////////////////////////////////////////////

TEST(Topaz2DReport_TemporalConvectionBC, RampedAmbient)
{
    SCOPED_TRACE("Begin Test: unit slab, linearly ramping ambient temperature.");

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

    constexpr auto dTime = 0.01;
    constexpr std::size_t nSteps = 200;
    constexpr auto rampRate = 1.0;
    constexpr auto hc = 1.0;

    std::vector<HygroThermFEM::FixedBCHCCoefficients> bcCoefficients;
    bcCoefficients.reserve(nSteps);
    for(std::size_t step = 1; step <= nSteps; ++step)
    {
        const double tAir = rampRate * static_cast<double>(step) * dTime;
        bcCoefficients.emplace_back(tAir, hc);
    }
    multiDomain.thermal().createBC_FixedHc(21, 22, bcCoefficients);

    auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
    auto humidities = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);
    std::vector<std::vector<double>> history;
    history.reserve(nSteps);
    for(std::size_t step = 0; step < nSteps; ++step)
    {
        const auto stepResult = multiDomain.transient(temperatures, humidities, dTime, step);
        temperatures = stepResult.temperature;
        humidities = stepResult.humidity;
        history.push_back(temperatures);
    }

    for(std::size_t step = 0; step < history.size(); ++step)
    {
    }

    // Exact series at x = 0, 0.5, 1.0 for steps 10, 50, 100, 200 (t = step * dt).
    const std::vector<std::pair<std::size_t, std::array<double, 3>>> checkpoints{
      {10, {0.000133270, 0.001741519, 0.019596747}},
      {50, {0.044252751, 0.074155792, 0.181104565}},
      {100, {0.221263680, 0.280553931, 0.470397249}},
      {200, {0.844065155, 0.937719514, 1.224394004}},
    };

    for(const auto & [step, expected] : checkpoints)
    {
        for(std::size_t pos = 0; pos < expected.size(); ++pos)
        {
            EXPECT_NEAR(expected[pos], history[step - 1][pos * 10], 0.005)
              << "step " << step << ", position " << pos;
        }
    }
}
