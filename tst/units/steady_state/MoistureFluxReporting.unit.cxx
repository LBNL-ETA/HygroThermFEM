#include <cmath>
#include <vector>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

#include "BeamBuilder.hxx"
#include "TestMaterials.hxx"

// The water flux the moisture domain reports must be the flux its equations transport. The
// vapour flux is -delta grad(c_sat(T) phi): a strip at uniform humidity under a temperature
// gradient carries a vapour flux through the gradient of c_sat alone, and the reported flux
// has to show it. Until 2026-11-03 the report differentiated phi only and gave zero here,
// and at a warm face of a real wall the wrong sign.
TEST(MoistureFluxReporting, UniformHumidityUnderATemperatureGradientCarriesVapor)
{
    HygroThermFEM::SimulationProperties::Instance().setCalculationParameters(
      true, true, true, true, false);
    HygroThermFEM::MultiDomain multiDomain;
    const auto & sandstone =
      multiDomain.materials().createSolidMaterial(TestHelper::CottaerSandstone());
    constexpr double thickness{0.1};
    constexpr std::size_t elements{20u};
    constexpr double humidity{0.5};
    TestHelper::BeamBuilder builder(multiDomain);
    builder.xStart(0.0)
      .height(thickness / static_cast<double>(elements))
      .numElementsY(1)
      .state({.temperature = 20.0, .humidity = humidity, .pressure = 101325.0, .liquidPercent = 1.0})
      .addSegment({.material = sandstone.name(), .numElementsX = elements, .width = thickness})
      .build();

    // A linear temperature profile, 20 C at the left end and 0 C at the right, on the nodes:
    // the builder numbers them column by column, two per column.
    const auto columns = builder.xCoordinates();
    auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
    for(std::size_t column = 0u; column < columns.size(); ++column)
    {
        const double temperature{20.0 - 200.0 * columns[column]};
        temperatures[2u * column] = temperature;
        temperatures[2u * column + 1u] = temperature;
    }
    multiDomain.nodes().updateNodeTemperatures(temperatures, true);

    const auto flux = multiDomain.moisture().flux();
    ASSERT_EQ(flux.size(), temperatures.size());

    // Expected at an interior column: the mean of the two adjacent elements' fluxes,
    // -delta(T) d(c_sat phi)/dx, with delta = D_air(T) / mu at the element's mean temperature.
    constexpr double resistanceFactor{15.0};
    const auto elementFlux = [&](const std::size_t left, const std::size_t right) {
        const double meanTemperature{(temperatures[2u * left] + temperatures[2u * right]) / 2.0};
        const double permeability{
          HygroThermFEM::vaporDiffusionCoefficientAtTemperature(meanTemperature) / resistanceFactor};
        const double concentrationLeft{
          HygroThermFEM::saturationConcentrationAtTemperature(temperatures[2u * left]) * humidity};
        const double concentrationRight{
          HygroThermFEM::saturationConcentrationAtTemperature(temperatures[2u * right]) * humidity};
        return -permeability * (concentrationRight - concentrationLeft) / (columns[right] - columns[left]);
    };
    for(std::size_t column = 1u; column + 1u < columns.size(); ++column)
    {
        const double expected{(elementFlux(column - 1u, column) + elementFlux(column, column + 1u)) / 2.0};
        EXPECT_GT(expected, 0.0) << "vapor moves toward the cold end, +x";
        EXPECT_NEAR(flux[2u * column].x, expected, 0.03 * expected) << "column " << column;
        EXPECT_NEAR(flux[2u * column].y, 0.0, 0.03 * expected) << "column " << column;
    }

    HygroThermFEM::SimulationProperties::Instance().resetCalculationParameters();
}
