#pragma once

#include <functional>
#include <string>
#include <vector>

#include "ErrorEstimator.hxx"

//! Regression fixtures captured from the legacy estimator, baked directly into
//! C++ (no external files, no legacy format). Each fixture lives in its own
//! tst/fixtures/<model>.cpp, builds an Input plus the legacy "golden" results,
//! and self-registers. The Regression test replays every registered case
//! through the new estimator and asserts the results match.
//!
//! Fixtures are produced once by converting a legacy `.eedump` (written by the
//! instrumented Common/errorest) into a baked .cpp; the ugly legacy format never
//! enters this module.
namespace lbnl::errorest::test
{
    struct LegacyCase
    {
        std::string name;
        Input input;
        double globalErrorPercent{0.0};      //!< legacy global error %
        std::vector<double> elementError{};  //!< legacy per-element relative error
        std::vector<int> refineElements{};   //!< element indices the legacy flagged
    };

    using Maker = std::function<LegacyCase()>;

    //! The set of registered fixtures (function-local static; safe across TUs).
    std::vector<Maker> & registry();

    //! Static-init self-registration: declare one per fixture .cpp.
    struct Registrar
    {
        explicit Registrar(Maker maker);
    };
}
