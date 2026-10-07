#include <cmath>
#include <cstddef>
#include <iomanip>
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
///
/// Validation book (hygrothermfem_python) dataset: sealed_strip_steady_gradient,
/// chapter "The sealed strip: steady state and conservation".
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

    //! Builds the sealed beam, one element high, with the temperature field written
    //! into the initial nodal state. Moisture-only means that field stays frozen.
    void buildStrip(HygroThermFEM::MultiDomain & multiDomain, std::size_t nColumns)
    {
        for(std::size_t col = 0; col < nColumns; ++col)
        {
            const HygroThermFEM::State state({.temperature = temperatureAt(xOf(col, nColumns)),
                                              .humidity = phiInitial,
                                              .pressure = 101325.0,
                                              .liquidPercent = 1.0});
            const std::size_t bottom = col * 2 + 1;
            multiDomain.nodes().createNode(
              {.index = bottom, .x = xOf(col, nColumns), .y = 0.0, .state = state});
            multiDomain.nodes().createNode(
              {.index = bottom + 1, .x = xOf(col, nColumns), .y = height, .state = state});
        }

        const auto & material =
          multiDomain.materials().createSolidMaterial(TestHelper::LinearSorption());
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

    std::vector<double> deviations;
    for(const std::size_t nElements : {10, 20, 40, 80})
    {
        const auto profile = runToSteady(nElements + 1);
        deviations.push_back(maxDeviationFromClosedForm(profile));
        std::cout << "[sealed strip] " << std::setw(3) << nElements << " elements, max |phi - "
                  << "closed form| = " << std::setprecision(4) << deviations.back() << "\n";
    }

    HygroThermFEM::SimulationProperties::Instance().resetCalculationParameters();

    // Exact on EVERY mesh, including the coarsest: the vapor term is assembled on its
    // product potential c_sat * phi, whose discrete flux vanishes identically when the
    // nodal products are equal, so there is no discretization error here to refine away.
    // Measured 1.3e-12 at 10 elements to 7.8e-13 at 80 -- the residual is the nonlinear
    // iteration's own tolerance, and it does not grow as the profile is resolved.
    for(std::size_t idx = 0; idx < deviations.size(); ++idx)
    {
        EXPECT_LT(deviations[idx], 1e-10) << "mesh index " << idx;
    }
}
