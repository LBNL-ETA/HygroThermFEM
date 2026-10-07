#include <cmath>
#include <cstddef>
#include <iomanip>
#include <vector>

#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "TestHelpers.hxx"
#include "TestMaterials.hxx"

/////////////////////////////////////////////////////////////////////////////////////
/// Conservation on a sealed strip, over ten simulated days, with the steep part of
/// the sorption isotherm in play.
///
/// Cottaer sandstone, 0.1 m, sealed on every boundary, humidity graded 0.6 -> 0.99
/// across the strip and an imposed 40 -> 20 C gradient driving vapour from the warm
/// end to the cold one. A closed system exchanges nothing, so the total stored water
/// is an invariant of the continuous problem, and the discrete operator inherits it
/// exactly: the transport matrices have vanishing column sums, so whatever leaves one
/// node arrives at another.
///
/// The measured drift is 1.05e-7 relative over the 240 steps. Its origin is NOT the
/// nonlinear iteration: sweeping the solver's error tolerance from 1e-6 to 1e-12
/// leaves the number bit-identical, as does removing the vapour storage term by
/// running the non-porous twin material. The 1D reference solver on the same case
/// stays at 4.7e-12. What sets the engine's level is not yet identified; the value
/// is pinned here so that it cannot move unnoticed.
///
/// Validation book (hygrothermfem_python) dataset: sealed_strip_conservation, chapter
/// "The sealed strip: steady state and conservation".
/////////////////////////////////////////////////////////////////////////////////////

namespace
{
    constexpr std::size_t nColumns = 21;   // 20 elements along x
    constexpr double length = 0.1;
    constexpr double height = 0.05;
    constexpr double phiLeft = 0.6;
    constexpr double phiRight = 0.99;
    constexpr double dTime = 3600.0;
    constexpr int nSteps = 240;

    double xOf(std::size_t col)
    {
        return length * static_cast<double>(col) / static_cast<double>(nColumns - 1);
    }

    double temperatureAt(double xPos)
    {
        return 40.0 - 20.0 * (xPos / length);
    }

    double humidityAt(std::size_t col)
    {
        const double fraction =
          static_cast<double>(col) / static_cast<double>(nColumns - 1);
        return phiLeft + (phiRight - phiLeft) * fraction;
    }

    void buildStrip(HygroThermFEM::MultiDomain & multiDomain)
    {
        for(std::size_t col = 0; col < nColumns; ++col)
        {
            const HygroThermFEM::State state({.temperature = temperatureAt(xOf(col)),
                                              .humidity = humidityAt(col),
                                              .pressure = 101325.0,
                                              .liquidPercent = 1.0});
            const std::size_t bottom = col * 2 + 1;
            multiDomain.nodes().createNode(
              {.index = bottom, .x = xOf(col), .y = 0.0, .state = state});
            multiDomain.nodes().createNode(
              {.index = bottom + 1, .x = xOf(col), .y = height, .state = state});
        }

        const auto & material =
          multiDomain.materials().createSolidMaterial(TestHelper::CottaerSandstone());
        for(std::size_t col = 0; col + 1 < nColumns; ++col)
        {
            const std::size_t left = col * 2 + 1;
            const std::size_t right = (col + 1) * 2 + 1;
            multiDomain.createElement({.node1 = left,
                                       .node2 = right,
                                       .node3 = right + 1,
                                       .node4 = left + 1,
                                       .material = material.name()});
        }
    }

    //! Total stored water of the beam, by the same half-cell assembly the reference
    //! solver uses: each element contributes its mean nodal water content.
    double totalWater(const HygroThermFEM::MultiDomain & multiDomain)
    {
        const auto water = multiDomain.nodes().properties(HygroThermFEM::Variable::water);
        const double dxElem = length / static_cast<double>(nColumns - 1);
        double total = 0.0;
        for(std::size_t col = 0; col + 1 < nColumns; ++col)
        {
            const std::size_t left = col * 2;
            const std::size_t right = (col + 1) * 2;
            const double mean = 0.25
                                * (water[left] + water[left + 1] + water[right]
                                   + water[right + 1]);
            total += mean * dxElem * height;
        }
        return total;
    }
}   // namespace

TEST(SealedStrip_Conservation, TotalWaterIsInvariant)
{
    SCOPED_TRACE("Begin Test: sealed strip conservation over ten days.");

    HygroThermFEM::MultiDomain multiDomain({.performThermal = false, .performMoisture = true});
    buildStrip(multiDomain);

    auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
    auto humidities = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);

    const double initialWater = totalWater(multiDomain);
    double maxDrift = 0.0;
    for(int step = 0; step < nSteps; ++step)
    {
        const auto solution = multiDomain.transient(temperatures, humidities, dTime, step);
        temperatures = solution.temperature;
        humidities = solution.humidity;
        maxDrift =
          std::max(maxDrift, std::abs(totalWater(multiDomain) - initialWater) / initialWater);
    }

    std::cout << "[sealed strip] relative drift in total water over " << nSteps
              << " steps = " << std::setprecision(4) << maxDrift << "\n";

    // Pinned at twice the measured 1.053e-7. The transport operator conserves
    // structurally (its column sums vanish), so this residue enters through the storage
    // side of the step rather than through the flux; see the header note.
    EXPECT_LT(maxDrift, 2.2e-7);
}
