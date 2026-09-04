#include <algorithm>
#include <cmath>
#include <vector>

#include <gtest/gtest.h>

#include "EN15026Material.hxx"
#include "HygroThermFEM2D.hxx"
#include "SlabCreator.hxx"
#include "TestHelpers.hxx"

/////////////////////////////////////////////////////////////////////////////////////
/// EN 15026:2007 Annex A -- the standard's liquid transport, validated on its own
/// against an exact similarity property rather than against stored numbers.
///
/// The annex specifies liquid transport as a liquid conductivity K(w), which the
/// engine consumes as a moisture diffusivity through the standard's eq. (19),
/// D_w = -K dp_suc/dw. That conversion is the part worth checking independently,
/// and it can be checked without any reference data at all.
///
/// Suppress vapour transport (a diffusion resistance factor of 1e9 makes the vapour
/// permeability negligible) and the model on a half line reduces to
///
///     dw/dt = d/dx( D_w(w) dw/dx )
///
/// with a constant initial state and a constant value imposed at x = 0. That problem
/// admits a similarity solution: w depends on x and t only through x / sqrt(t), for
/// ANY shape of D_w(w). So the profile at time t and the profile at 4t must lie on top
/// of each other once the abscissa is doubled, and they must do so pointwise.
///
/// This catches the two ways the conversion can be wrong, separately. A D_w of the
/// wrong SHAPE breaks the collapse, because the profiles are then not related by a
/// pure stretch. A D_w of the wrong SCALE keeps the shape but moves the front at the
/// wrong rate, which the same comparison sees as a horizontal offset. Neither needs
/// the standard's tables, so this test stands even though the engine cannot yet run
/// the full benchmark (its material model carries one scalar diffusion resistance
/// factor while the annex's is moisture dependent -- see EN15026_UptakeProjection).
///
/// Mirrors the reference solver's test_liquid_only_uptake_collapses_in_boltzmann_variable
/// (hygrothermfem_python), on the same material, mesh rule and schedule.
/////////////////////////////////////////////////////////////////////////////////////

namespace
{
    //! Linear interpolation of a bottom-row profile at position x.
    double profileAt(const std::vector<double> & coords,
                     const std::vector<double> & values,
                     const double position)
    {
        const auto upper = std::ranges::upper_bound(coords, position) - coords.begin();
        const auto right = (std::min)(static_cast<std::size_t>(upper), coords.size() - 1);
        const auto left = right - 1;
        const double fraction = (position - coords[left]) / (coords[right] - coords[left]);
        return values[left] + fraction * (values[right] - values[left]);
    }
}   // namespace

TEST(EN15026_LiquidUptakeSimilarity, ProfilesCollapseInTheBoltzmannVariable)
{
    SCOPED_TRACE("Begin Test: EN 15026 liquid-only uptake, similarity in x / sqrt(t).");

    constexpr double initialHumidity = 0.5;
    constexpr double surfaceHumidity = 0.95;
    constexpr double temperature = 20.0;
    constexpr double dTime = 5400.0;
    constexpr std::size_t earlyStep = 256u;
    constexpr std::size_t lateStep = 1024u;   //!< four times earlyStep, so the abscissa doubles

    // Vapour transport suppressed, so only the annex's liquid curve moves moisture.
    HygroThermFEM::MultiDomain multiDomain({.performThermal = false, .performMoisture = true});
    const auto & material =
      multiDomain.materials().createSolidMaterial(TestHelper::EN15026AnnexA(1.0e9));

    const auto coords = TestHelper::en15026GradedCoordinates(1.0e-4, 1.04, 0.5);
    const auto nColumns = coords.size();

    TestHelper::SlabBuilder(multiDomain)
      .gridXCoordinates(coords)
      .height(0.05)
      .material(material.name())
      .state({.temperature = temperature,
              .humidity = initialHumidity,
              .pressure = 101325.0,
              .liquidPercent = 1.0})
      .build();

    multiDomain.moisture().createBC_FixedHumidity(
      1, 2, HygroThermFEM::TemperatureAndHumidity{temperature, surfaceHumidity});

    auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
    auto humidities = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);
    const auto initialWater =
      TestHelper::bottomRow(multiDomain.nodes().properties(HygroThermFEM::Variable::water),
                            nColumns,
                            2);

    std::vector<double> early;
    std::vector<double> late;
    for(std::size_t step = 1u; step <= lateStep; ++step)
    {
        const auto solution = multiDomain.transient(temperatures, humidities, dTime, step - 1u);
        temperatures = solution.temperature;
        humidities = solution.humidity;
        if(step == earlyStep)
        {
            early = TestHelper::bottomRow(solution.waterContent, nColumns, 2);
        }
        if(step == lateStep)
        {
            late = TestHelper::bottomRow(solution.waterContent, nColumns, 2);
        }
    }
    ASSERT_EQ(nColumns, early.size());
    ASSERT_EQ(nColumns, late.size());

    // The front has genuinely moved into the material, so the collapse below is a
    // statement about a real profile rather than about a flat field.
    EXPECT_GT(profileAt(coords, early, 0.005), 100.0) << "the wetting front did not advance";

    // The far end is still at its initial water content, so the half-line assumption the
    // similarity solution rests on still holds at the later time.
    EXPECT_NEAR(late.back(), initialWater.back(), 1e-9) << "the far end was disturbed";

    // The similarity itself: w(x, 4t) against w(2x, t). Measured worst difference 0.25 kg/m3
    // here, against 0.11 in the 1D reference solver, over a profile spanning about 86 kg/m3.
    // So the collapse holds to about a third of a percent of the range; the engine is looser
    // than the reference by roughly a factor of two, which is the discretisation difference
    // between a 2D quadrilateral mesh with lumped capacity and a 1D solver, not a property
    // error. A conversion mistake would miss by a large fraction of the range, not by this.
    constexpr double collapseTolerance = 0.5;
    double worstCollapse{0.0};
    double worstPosition{0.0};
    for(std::size_t sample = 0u; sample < 200u; ++sample)
    {
        const double position = 0.04 * static_cast<double>(sample) / 199.0;
        const double difference =
          std::abs(profileAt(coords, late, 2.0 * position) - profileAt(coords, early, position));
        if(difference > worstCollapse)
        {
            worstCollapse = difference;
            worstPosition = position;
        }
    }
    EXPECT_LT(worstCollapse, collapseTolerance)
      << "profiles at t and 4t do not collapse in x / sqrt(t): worst difference "
      << worstCollapse << " kg/m3 at x = " << worstPosition << " m, which means the annex's "
      << "K to D_w conversion has the wrong shape or the wrong scale";
}
