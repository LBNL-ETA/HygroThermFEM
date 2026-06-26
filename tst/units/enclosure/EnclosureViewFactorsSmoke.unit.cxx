#include <vector>
#include <gtest/gtest.h>

#include <WCEViewer.hpp>

// B0 smoke test: confirm HygroThermFEM builds against the WCE feature branch (via the
// KeffCavity -> Windows-CalcEngine local-sibling cascade) and can call the new enclosure
// view-factor engine end-to-end. A unit square gives the known closed-enclosure view factors.
TEST(EnclosureViewFactorsSmoke, UnitSquareThroughWCE)
{
    using namespace Viewer;

    const std::vector<RadiationSegment> segments{
      {.startPoint = {0.0, 0.0}, .endPoint = {0.0, 1.0}, .emissivity = 0.9, .enclosureId = 0u},
      {.startPoint = {0.0, 1.0}, .endPoint = {1.0, 1.0}, .emissivity = 0.9, .enclosureId = 0u},
      {.startPoint = {1.0, 1.0}, .endPoint = {1.0, 0.0}, .emissivity = 0.9, .enclosureId = 0u},
      {.startPoint = {1.0, 0.0}, .endPoint = {0.0, 0.0}, .emissivity = 0.9, .enclosureId = 0u}};

    const auto result = computeEnclosureViewFactors(segments, {});

    EXPECT_EQ(4u, result.viewFactors.size());
    EXPECT_NEAR(0.292893219, result.viewFactors(0, 1), 1e-6);
    EXPECT_NEAR(0.414213562, result.viewFactors(0, 2), 1e-6);
    EXPECT_NEAR(0.292893219, result.viewFactors(0, 3), 1e-6);

    // Closed enclosure: each segment's environment (deficit) view factor is ~0.
    for(const auto environment : result.environmentViewFactor)
    {
        EXPECT_NEAR(0.0, environment, 1e-9);
    }
}
