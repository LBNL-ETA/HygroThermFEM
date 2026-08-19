#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <vector>

#include <gtest/gtest.h>

#include "ErrorEstimator.hxx"

#include "RegressionFixture.hxx"

#ifndef ERROREST_FIXTURE_DIR
#define ERROREST_FIXTURE_DIR ""
#endif

using namespace HygroThermFEM::errorest;

namespace fs = std::filesystem;

namespace
{
    double tol(double magnitude)
    {
        // The new SPR solve and legacy esterror solve the same patch normal equations, so
        // they agree to ~1e-15 on well-conditioned patches. On near-rank-deficient sliver
        // patches (nearly collinear points) both carry inherent ill-conditioning round-off
        // that differs by up to ~1e-4 between the two elimination paths -- benign: the
        // global error and the refine set still match exactly. 1e-4 covers it while staying
        // ~100x below any algorithmic difference (material/region mis-grouping is ~1e-2).
        return 1e-4 * (1.0 + std::abs(magnitude));
    }

    std::vector<fs::path> fixtureFiles()
    {
        std::vector<fs::path> files;
        const fs::path dir{ERROREST_FIXTURE_DIR};
        if (fs::exists(dir))
        {
            for (const auto & entry : fs::directory_iterator(dir))
            {
                if (entry.is_regular_file() && entry.path().extension() == ".json")
                {
                    files.push_back(entry.path());
                }
            }
        }
        std::ranges::sort(files);
        return files;
    }
}

//! Replays every captured legacy fixture (JSON in tst/units/errorest/fixtures/)
//! through the new estimator and asserts the results match. Add a case by dropping
//! a <model>.json fixture into that directory — no recompile of this file needed.
TEST(Regression, MatchesLegacy)
{
    const auto files = fixtureFiles();
    if (files.empty())
    {
        GTEST_SKIP() << "no .json fixtures in " << ERROREST_FIXTURE_DIR;
    }

    for (const auto & file : files)
    {
        const auto name = file.filename().string();
        const auto golden = test::loadFixture(file.string());
        const auto result = estimate(golden.input);
        ASSERT_TRUE(result.has_value()) << name;

        EXPECT_NEAR(result->globalErrorPercent, golden.globalErrorPercent,
                    tol(golden.globalErrorPercent))
            << name << ": global error";

        ASSERT_EQ(result->elementRelativeError.size(), golden.elementError.size()) << name;
        for (std::size_t idx = 0; idx < golden.elementError.size(); ++idx)
        {
            EXPECT_NEAR(result->elementRelativeError[idx], golden.elementError[idx],
                        tol(golden.elementError[idx]))
                << name << ": element " << idx;
        }

        auto refine = result->elementsToRefine;
        std::ranges::sort(refine);
        EXPECT_EQ(refine, golden.refineElements) << name << ": refinement set";
    }
}
