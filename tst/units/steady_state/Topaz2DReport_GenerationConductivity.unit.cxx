#include <gtest/gtest.h>

#include "DumpCsv.hxx"
#include "HygroThermFEM2D.hxx"
#include "SlabCreator.hxx"
#include "TestMaterials.hxx"

/////////////////////////////////////////////////////////////////////////////////////
/// TOPAZ2D Validation Status Report (Drayton, UCRL-ID-106550, December 1990),
/// section 2.9a -- Plane Steady Nonlinear, Material Type 3 (k(T)).
///
/// Unit slab with uniform generation q = 1 AND linear temperature-dependent
/// conductivity k(T) = k_a (1 + beta T), beta = 1, adiabatic at x = 0, surface at
/// x = 1 held at 0. Closed form (Arpaci p. 130):
///
///     (T - Ts)/S = (-1 + sqrt(1 + 2 beta S [1 - (x/L)^2])) / (beta S),
///     S = q L^2 / (2 k_s)
///
/// so T(0) = sqrt(2) - 1 = 0.414214 (the report's published table). The k(T)
/// dependence enters through the tabular temperature-dependent conductivity --
/// exact for a linear law, since piecewise-linear interpolation of a linear
/// function is that function -- enabled via the global simulation flag, which the
/// element reads at construction time. The steady coupling loop re-linearizes
/// k about each pass's temperature (the same Picard treatment as radiation).
///
/// This is the report's only feasible temperature-dependent-conductivity problem:
/// its closed form requires the generation term.
/////////////////////////////////////////////////////////////////////////////////////

TEST(Topaz2DReport_GenerationConductivity, MaterialType3)
{
    SCOPED_TRACE("Begin Test: steady slab, uniform generation, k(T) = 1 + T.");

    HygroThermFEM::SimulationProperties::Instance().setCalculationParameters(
      false, false, false, false, true);
    HygroThermFEM::SimulationProperties::Instance().setIterationParameters(1.0, 1e-10, 500);

    {
        HygroThermFEM::MultiDomain multiDomain(
          {.performThermal = true, .performMoisture = false});

        constexpr HygroThermFEM::State state({
            .temperature = 0.0,
            .humidity = 0.0,
            .pressure = 101325.0,
            .liquidPercent = 0
        });

        auto params = TestHelper::TestMaterial();
        params.density = 1.0;
        params.heatCapacity = 1.0;
        // k(T) = 1 + T over the solution range [0, 0.42]; the moisture curve is flat.
        params.thermalConductivityTemperatureDependent = {{0.0, 1.0}, {1.0, 2.0}};
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

        TestHelper::CsvDump dump("topaz_generation_kt.csv", 11);
        dump.addRow(1, TestHelper::bottomRow(solution.temperature, 11, 2));

        // The report's published analytic column (p. 17), T(0) = sqrt(2) - 1.
        const std::vector<double> expected{0.414214,
                                           0.410674,
                                           0.400000,
                                           0.382027,
                                           0.356466,
                                           0.322876,
                                           0.280625,
                                           0.228821,
                                           0.166190,
                                           0.0908712,
                                           0.0};

        for(std::size_t col = 0; col < expected.size(); ++col)
        {
            EXPECT_NEAR(expected[col], solution.temperature[col * 2], 1e-4)
              << "column " << col;
        }
    }

    HygroThermFEM::SimulationProperties::Instance().resetIterationParameters();
    HygroThermFEM::SimulationProperties::Instance().resetCalculationParameters();
}
