#include <algorithm>
#include <cmath>
#include <cstddef>

#include <gtest/gtest.h>

#include "ErrorEstimator.hxx"

#include "RegressionFixture.hxx"

using namespace lbnl::errorest;

namespace
{
    double tol(double magnitude)
    {
        // The new SPR solve (WCE LUFactor, no pivoting) and the legacy ludcmp
        // (partial pivoting) agree up to floating-point rounding for the SPD
        // patch normal equations, so a tight relative tolerance is expected.
        return 1e-6 * (1.0 + std::abs(magnitude));
    }
}

//! Replays every baked legacy fixture through the new estimator and asserts the
//! results match. Add a case by dropping a fixture .cpp into tst/fixtures/.
TEST(Regression, MatchesLegacy)
{
    const auto & makers = test::registry();
    if (makers.empty())
    {
        GTEST_SKIP() << "no baked legacy fixtures yet (convert a .eedump into "
                        "tst/fixtures/<model>.cpp to enable this test)";
    }

    for (const auto & make : makers)
    {
        const auto golden = make();
        const auto result = estimate(golden.input);
        ASSERT_TRUE(result.has_value()) << golden.name;

        EXPECT_NEAR(result->globalErrorPercent, golden.globalErrorPercent,
                    tol(golden.globalErrorPercent))
            << golden.name << ": global error";

        ASSERT_EQ(result->elementRelativeError.size(), golden.elementError.size()) << golden.name;
        for (std::size_t idx = 0; idx < golden.elementError.size(); ++idx)
        {
            EXPECT_NEAR(result->elementRelativeError[idx], golden.elementError[idx],
                        tol(golden.elementError[idx]))
                << golden.name << ": element " << idx;
        }

        auto refine = result->elementsToRefine;
        std::ranges::sort(refine);
        EXPECT_EQ(refine, golden.refineElements) << golden.name << ": refinement set";
    }
}
