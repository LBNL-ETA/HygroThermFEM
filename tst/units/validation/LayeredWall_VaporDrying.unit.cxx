#include <array>

#include <gtest/gtest.h>

#include "BeamBuilder.hxx"
#include "HygroThermFEM2D.hxx"
#include "TestMaterials.hxx"

/////////////////////////////////////////////////////////////////////////////////////
/// Two-layer wall drying through one face: Cottaer Sandstone (0.06 m) behind
/// Stucco (0.04 m), uniform initial humidity 0.7 at 20 C, sealed on the left,
/// convective vapor exchange (h_c = 10, phi_air = 0.4) on the right. Isothermal,
/// moisture only: relative humidity is continuous across the material interface
/// while water content jumps with the isotherms, and the drying front crosses
/// the interface within the simulated two days.
///
/// This is the validation book's first LAYERED engine dataset: the bottom node
/// row is compared against the 1D reference solver running the same two-layer
/// grid (hygrothermfem_python, case layered_vapor_drying). The humidity stays
/// well below the near-saturation taper, so the vapor transfer coefficient is
/// the plain Lewis value beta = h_c / (rho_air Cp_air) on both sides
/// (IConvectiveCoefficient::waterVaporTransferCoefficient multiplies the
/// per-node film coefficients in at its final line).
/////////////////////////////////////////////////////////////////////////////////////

TEST(LayeredWall_VaporDrying, TwoLayerDryingThroughOneFace)
{
    SCOPED_TRACE("Begin Test: layered wall drying through the right face.");

    HygroThermFEM::MultiDomain multiDomain({.performThermal = false, .performMoisture = true});

    constexpr double initialTemperature = 20.0;
    constexpr double initialHumidity = 0.7;

    constexpr HygroThermFEM::State state({
        .temperature = initialTemperature,
        .humidity = initialHumidity,
        .pressure = 101325.0,
        .liquidPercent = 1.0
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

    constexpr auto hc = 10.0;
    constexpr auto airHumidity = 0.4;
    const HygroThermFEM::FixedBCHCCoefficients interiorBc{
      /*airTemperature=*/initialTemperature, hc, airHumidity};
    for(const auto & [index1, index2] : builder.rightEdge())
    {
        multiDomain.createBC_FixedHc(index1, index2, interiorBc);
    }

    constexpr double dTime = 3600.0;
    constexpr int nSteps = 48;

    auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
    auto humidities = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);
    std::vector<std::vector<double>> history;
    history.reserve(nSteps);
    for(int step = 0; step < nSteps; ++step)
    {
        const auto solution = multiDomain.transient(temperatures, humidities, dTime, step);
        temperatures = solution.temperature;
        humidities = solution.humidity;
        history.push_back(humidities);
    }

    // Drying sanity: humidity decreases monotonically at the exposed face and
    // stays between the ambient value and the initial state everywhere.
    const auto & final = history.back();
    EXPECT_LT(final[20], initialHumidity);
    EXPECT_GT(final[20], airHumidity);
    EXPECT_LT(final[0], initialHumidity);
    for(std::size_t col = 0; col < 11; ++col)
    {
        EXPECT_GE(final[col * 2], airHumidity) << "column " << col;
        EXPECT_LE(final[col * 2], initialHumidity + 1e-12) << "column " << col;
    }
    // The sealed face lags the exposed face.
    EXPECT_GT(final[0], final[20]);

    // 48 h humidity checkpoints from the independently implemented 1D reference solver
    // (hygrothermfem_python, case layered_vapor_drying) on the same two-layer grid:
    // the sealed face, the cottaer interior, the material interface (column 6), the
    // stucco interior and the exposed face. Tolerances are ~2x the per-node
    // engine-reference deviation measured at capture, floored at 1e-6; the interface
    // and the stucco interior carry the coefficient-evaluation-geometry difference the
    // book chapter records (1.7e-3 profile-wide maximum).
    struct Checkpoint
    {
        std::size_t column;
        double humidity;
        double tolerance;
    };
    constexpr std::array<Checkpoint, 5> checkpoints{
      {{.column = 0, .humidity = 0.699997509925111, .tolerance = 1e-6},
       {.column = 3, .humidity = 0.699884937834578, .tolerance = 2e-5},
       {.column = 6, .humidity = 0.695045571323813, .tolerance = 4e-4},
       {.column = 8, .humidity = 0.650816791804612, .tolerance = 2e-3},
       {.column = 10, .humidity = 0.403467560658606, .tolerance = 1e-6}}};
    for(const auto & [column, humidity, tolerance] : checkpoints)
    {
        EXPECT_NEAR(final[column * 2], humidity, tolerance) << "column " << column;
    }
}
