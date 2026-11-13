#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "TestHelpers.hxx"
#include "TestMaterials.hxx"

/////////////////////////////////////////////////////////////////////////////////////
/// Sealed strip at steady state under an imposed 40 -> 20 C gradient, checked
/// against the closed form.
///
/// A 0.1 m strip with a linear isotherm (w = phi), sealed on every boundary,
/// starts at a uniform humidity and is held in a stationary linear temperature
/// field. Vapor transport alone redistributes moisture until the flux vanishes.
///
/// At that steady state the flux is zero EVERYWHERE, not merely its divergence,
/// so phi * c_sat(T) is constant along the strip:
///
///     phi(x) = C / c_sat(T(x))
///
/// and C follows from conservation, the strip being closed. The expected profile
/// is therefore evaluated here from the closed form rather than pasted in as
/// constants, using the engine's own saturation function for c_sat.
///
/// The steady state is reached to machine precision (the profile stops changing
/// between steps), and total water is conserved to 1e-13, so what the deviation
/// from the closed form measures is the discretization of the vapor
/// temperature-gradient term alone: it falls at second order under refinement.
/////////////////////////////////////////////////////////////////////////////////////

namespace
{
    constexpr double length = 0.1;
    constexpr double height = 0.05;
    constexpr double phiInitial = 0.4;
    constexpr double dTime = 36000.0;
    constexpr int nSteps = 300;

    double xOf(std::size_t col, std::size_t nColumns)
    {
        return length * static_cast<double>(col) / static_cast<double>(nColumns - 1);
    }

    //! Imposed field: 40 C at x = 0 falling linearly to 20 C at x = L.
    double temperatureAt(double xPos)
    {
        return 40.0 - 20.0 * (xPos / length);
    }

    //! Nodal control volumes of the 1D column layout: half an element on each side of
    //! a node, so the end columns carry half the interior value. The beam height is a
    //! common factor and cancels in the conservation ratio below.
    std::vector<double> lumpedVolumes(std::size_t nColumns)
    {
        const double dxElem = length / static_cast<double>(nColumns - 1);
        std::vector<double> volumes(nColumns, dxElem);
        volumes.front() = 0.5 * dxElem;
        volumes.back() = 0.5 * dxElem;
        return volumes;
    }

    //! Closed-form steady humidity profile, phi = C / c_sat(T), with C fixed by the
    //! initial water content of the sealed strip. With a linear isotherm w = phi the
    //! stored water is proportional to phi, so C = phi0 * sum(V) / sum(V / c_sat).
    std::vector<double> zeroFluxSteadyState(std::size_t nColumns)
    {
        const auto volumes = lumpedVolumes(nColumns);

        std::vector<double> saturation;
        saturation.reserve(nColumns);
        for(std::size_t col = 0; col < nColumns; ++col)
        {
            saturation.push_back(HygroThermFEM::saturationConcentrationAtTemperature(
              temperatureAt(xOf(col, nColumns))));
        }

        double totalVolume = 0.0;
        double weightedVolume = 0.0;
        for(std::size_t col = 0; col < nColumns; ++col)
        {
            totalVolume += volumes[col];
            weightedVolume += volumes[col] / saturation[col];
        }
        const double constant = phiInitial * totalVolume / weightedVolume;

        std::vector<double> exact;
        exact.reserve(nColumns);
        for(const double csat : saturation)
        {
            exact.push_back(constant / csat);
        }
        return exact;
    }

    //! The nodes of a strip, one element high, with the temperature field written into the
    //! initial nodal state. Moisture-only means that field stays frozen. The node indices
    //! start after firstIndex, so two strips fit in one model.
    void createStripNodes(HygroThermFEM::MultiDomain & multiDomain,
                          const std::size_t nColumns,
                          const std::size_t firstIndex = 0u)
    {
        for(std::size_t col = 0; col < nColumns; ++col)
        {
            const HygroThermFEM::State state({.temperature = temperatureAt(xOf(col, nColumns)),
                                              .humidity = phiInitial,
                                              .pressure = 101325.0,
                                              .liquidPercent = 1.0});
            const std::size_t bottom = firstIndex + col * 2 + 1;
            multiDomain.nodes().createNode(
              {.index = bottom, .x = xOf(col, nColumns), .y = 0.0, .state = state});
            multiDomain.nodes().createNode(
              {.index = bottom + 1, .x = xOf(col, nColumns), .y = height, .state = state});
        }
    }

    //! The elements of a strip over nodes already created. Every node of the model must
    //! exist before the first element does.
    void createStripElements(HygroThermFEM::MultiDomain & multiDomain,
                             const std::size_t nColumns,
                             const std::string & materialName,
                             const std::size_t firstIndex = 0u)
    {
        for(std::size_t col = 0; col + 1 < nColumns; ++col)
        {
            const std::size_t left = firstIndex + col * 2 + 1;
            const std::size_t right = firstIndex + (col + 1) * 2 + 1;
            multiDomain.createElement({.node1 = left,
                                       .node2 = right,
                                       .node3 = right + 1,
                                       .node4 = left + 1,
                                       .material = materialName});
        }
    }

    //! The strip's material, w = phi, added to the model; returns its name.
    std::string addStripMaterial(HygroThermFEM::MultiDomain & multiDomain)
    {
        return multiDomain.materials().createSolidMaterial(TestHelper::LinearSorption()).name();
    }

    //! Builds the sealed strip of the closed form: its material, nodes and elements.
    void buildStrip(HygroThermFEM::MultiDomain & multiDomain, const std::size_t nColumns)
    {
        createStripNodes(multiDomain, nColumns);
        createStripElements(multiDomain, nColumns, addStripMaterial(multiDomain));
    }

    //! Water of a strip profile over the water of the uniform start, minus one. The
    //! isotherm is w = phi, so the water is the humidity summed over the nodal volumes.
    double relativeWaterChange(const std::vector<double> & profile)
    {
        const auto volumes = lumpedVolumes(profile.size());
        double water = 0.0;
        double startWater = 0.0;
        for(std::size_t col = 0; col < profile.size(); ++col)
        {
            water += volumes[col] * profile[col];
            startWater += volumes[col] * phiInitial;
        }
        return water / startWater - 1.0;
    }

    //! Runs the sealed strip to steady state and returns the bottom node row.
    std::vector<double> runToSteady(std::size_t nColumns)
    {
        HygroThermFEM::MultiDomain multiDomain(
          {.performThermal = false, .performMoisture = true});
        buildStrip(multiDomain, nColumns);

        auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
        auto humidities = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);
        for(int step = 0; step < nSteps; ++step)
        {
            const auto solution = multiDomain.transient(temperatures, humidities, dTime, step);
            temperatures = solution.temperature;
            humidities = solution.humidity;
        }
        return TestHelper::bottomRow(humidities, nColumns, 2);
    }

    double maxDeviationFromClosedForm(const std::vector<double> & profile)
    {
        const auto exact = zeroFluxSteadyState(profile.size());
        double deviation = 0.0;
        for(std::size_t col = 0; col < profile.size(); ++col)
        {
            deviation = std::max(deviation, std::abs(profile[col] - exact[col]));
        }
        return deviation;
    }
}   // namespace

TEST(SealedStrip_SteadyGradient, ZeroFluxClosedForm)
{
    SCOPED_TRACE("Begin Test: sealed strip driven to its zero-flux steady state.");

    // Vapor transport alone: the closed form describes the vapor equation, and the
    // material's liquid curve is zero in any case.
    HygroThermFEM::SimulationProperties::Instance().setCalculationParameters(
      true, false, false, false, false);

    constexpr std::array<std::size_t, 4> meshSizes{10u, 20u, 40u, 80u};
    std::vector<double> deviations;
    for(const std::size_t nElements : meshSizes)
    {
        deviations.push_back(maxDeviationFromClosedForm(runToSteady(nElements + 1)));
    }

    HygroThermFEM::SimulationProperties::Instance().resetCalculationParameters();

    // Exact on EVERY mesh, including the coarsest: the vapor term is assembled on its
    // product potential c_sat * phi, whose discrete flux vanishes identically when the
    // nodal products are equal, so there is no discretization error here to refine away.
    // What remains is the nonlinear iteration's own tolerance, and it does not grow as the
    // profile is resolved: 1.3e-12 to 7.8e-13 with the lagged-coefficient solve, which was
    // exact for this linear case; 5e-11 to 1.3e-9 since the moisture Newton takes its
    // direction from a finite-difference Jacobian (2026-09-03), whose truncation error sets
    // the floor.
    for(std::size_t idx = 0; idx < deviations.size(); ++idx)
    {
        EXPECT_LT(deviations[idx], 1e-8) << meshSizes[idx] << " elements";
    }
}

//! The steady solve on the same strip. Nothing enters or leaves it, so its flux equations
//! fix the profile only up to the water; the solve keeps the water of the start field, and
//! with it lands on the closed form the transient run settles to.
TEST(SealedStrip_SteadyGradient, TheSteadySolveKeepsTheStartWater)
{
    SCOPED_TRACE("Begin Test: sealed strip solved as a steady state.");

    HygroThermFEM::SimulationProperties::Instance().setCalculationParameters(
      true, false, false, false, false);

    constexpr std::array<std::size_t, 4> meshSizes{10u, 20u, 40u, 80u};
    for(const std::size_t nElements : meshSizes)
    {
        HygroThermFEM::MultiDomain multiDomain({.performThermal = false, .performMoisture = true});
        buildStrip(multiDomain, nElements + 1);

        const auto solution = multiDomain.steadyState();
        const auto profile = TestHelper::bottomRow(solution.humidity, nElements + 1, 2);

        EXPECT_TRUE(solution.converged) << nElements << " elements";
        EXPECT_LT(maxDeviationFromClosedForm(profile), 1e-8) << nElements << " elements";
        EXPECT_LT(std::abs(relativeWaterChange(profile)), 1e-12) << nElements << " elements";
    }

    HygroThermFEM::SimulationProperties::Instance().resetCalculationParameters();
}

//! Two strips in one model, not connected. One is sealed and keeps its water; the other
//! has a film at each end and goes to the humidity of the air there. The condition that
//! holds the water belongs to the sealed strip alone.
TEST(SealedStrip_SteadyGradient, ASealedStripBesideAnExchangingOne)
{
    SCOPED_TRACE("Begin Test: sealed strip beside a strip with films.");

    HygroThermFEM::SimulationProperties::Instance().setCalculationParameters(
      true, false, false, false, false);

    constexpr std::size_t nColumns{21u};
    constexpr std::size_t secondStrip{2u * nColumns};
    constexpr double airHumidity{0.8};
    HygroThermFEM::MultiDomain multiDomain({.performThermal = false, .performMoisture = true});
    const auto materialName = addStripMaterial(multiDomain);
    createStripNodes(multiDomain, nColumns);
    createStripNodes(multiDomain, nColumns, secondStrip);
    createStripElements(multiDomain, nColumns, materialName);
    createStripElements(multiDomain, nColumns, materialName, secondStrip);
    // Each film sees air at its own end's temperature, so the air's vapour pressure differs
    // between the ends and the exchanging strip carries a flux: a state of its own.
    multiDomain.createBC_FixedHc(secondStrip + 1, secondStrip + 2, {40.0, 5.0, airHumidity});
    multiDomain.createBC_FixedHc(
      secondStrip + 2 * nColumns, secondStrip + 2 * nColumns - 1, {20.0, 5.0, airHumidity});

    const auto solution = multiDomain.steadyState();
    const std::vector<double> sealed(solution.humidity.begin(), solution.humidity.begin() + secondStrip);
    const std::vector<double> exchanging(solution.humidity.begin() + secondStrip, solution.humidity.end());
    const auto sealedProfile = TestHelper::bottomRow(sealed, nColumns, 2);
    const auto exchangingProfile = TestHelper::bottomRow(exchanging, nColumns, 2);

    HygroThermFEM::SimulationProperties::Instance().resetCalculationParameters();

    EXPECT_TRUE(solution.converged);
    EXPECT_LT(maxDeviationFromClosedForm(sealedProfile), 1e-8);
    EXPECT_LT(std::abs(relativeWaterChange(sealedProfile)), 1e-12);
    // The films pull both ends of the other strip towards the air's humidity, far from the
    // 0.4 it started with: its water is not kept.
    EXPECT_NEAR(exchangingProfile.front(), airHumidity, 0.05);
    EXPECT_NEAR(exchangingProfile.back(), airHumidity, 0.05);
    EXPECT_GT(relativeWaterChange(exchangingProfile), 0.5);
}
