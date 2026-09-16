#pragma once

#include <string>
#include <vector>

namespace HygroThermFEM
{
    class IMaterial;
    class Materials;

    //! \brief A moisture table keyed by water content that does not end where the
    //! material's moisture storage function does.
    //!
    //! Every such table is read at the water content the sorption curve returns, so its
    //! axis is expected to reach the curve's maximum. A shorter table is not an error --
    //! the solver extrapolates its last segment -- but the material is then evaluated on
    //! data it does not carry. Deliberately OUTSIDE the calculation path: no material or
    //! domain code calls this; a consumer (the mediator, a GUI) runs it when it wants to
    //! tell the user.
    struct WaterContentAxisMismatch
    {
        std::string materialName;
        std::string tableName;
        double tableMaximum{0.0};
        double sorptionMaximum{0.0};

        //! One user-facing sentence naming the material, the table and both maxima.
        [[nodiscard]] std::string message() const;
    };

    //! \brief Checks the moisture-dependent thermal conductivity, the liquid transport
    //! coefficient and the moisture-dependent vapor resistance factor against the
    //! material's sorption curve.
    //!
    //! A material without a sorption curve has nothing to compare against and yields no
    //! mismatch. A constant table -- one point, or points at a single water content -- holds
    //! at every water content and is never reported.
    [[nodiscard]] std::vector<WaterContentAxisMismatch>
      checkWaterContentAxes(const IMaterial & material);

    //! \brief The same check over every solid material in the pool, in pool order.
    [[nodiscard]] std::vector<WaterContentAxisMismatch>
      checkWaterContentAxes(const Materials & materials);
}   // namespace HygroThermFEM
