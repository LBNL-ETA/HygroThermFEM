#pragma once

#include <array>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <fileParse/Common.hxx>
#include <fileParse/Vector.hxx>
#include <fileParse/FileDataHandler.hxx>

#include "ErrorEstimator.hxx"

//! Loads an error-estimator regression fixture (one captured legacy case) and builds an
//! Input plus the legacy "golden" results. Serialization goes through FileParse (the house
//! library), so a fixture can be JSON or XML — Common::loadFromFile auto-detects the format.
//! FileParse is a test-only dependency; nothing here reaches the HygroThermFEM library.
//! See tst/units/errorest/fixtures/*.json.
namespace lbnl::errorest::test
{
    //! On-disk shape (mirrors the JSON/XML). Coordinates/fluxes are flat x,y pairs.
    struct FixtureElement
    {
        double conductivity{0.0};
        int region{-1};                     //!< Mesh region id; -1 if the fixture omits it.
        std::vector<double> vertices{};     //!< x0,y0,x1,y1,...
        std::vector<double> gaussPoints{};  //!< x0,y0,...,x3,y3
        std::vector<double> flux{};         //!< fx0,fy0,...,fx3,fy3
        double legacyError{0.0};
        bool refine{false};
    };

    struct FixtureCase
    {
        std::string name{};
        double errorLimitPercent{0.0};
        double globalErrorPercent{0.0};
        std::vector<FixtureElement> elements{};
    };

    template<typename NodeAdapter>
    inline const NodeAdapter & operator>>(const NodeAdapter & node, FixtureElement & ele)
    {
        using FileParse::Child;
        using FileParse::operator>>;

        node >> Child{"Conductivity", ele.conductivity};
        node >> Child{"Region", ele.region};
        node >> Child{{"Vertices", "Coord"}, ele.vertices};
        node >> Child{{"GaussPoints", "Coord"}, ele.gaussPoints};
        node >> Child{{"Flux", "Comp"}, ele.flux};
        node >> Child{"LegacyError", ele.legacyError};
        node >> Child{"Refine", ele.refine};

        return node;
    }

    template<typename NodeAdapter>
    inline const NodeAdapter & operator>>(const NodeAdapter & node, FixtureCase & fix)
    {
        using FileParse::Child;
        using FileParse::operator>>;

        node >> Child{"Name", fix.name};
        node >> Child{"ErrorLimitPercent", fix.errorLimitPercent};
        node >> Child{"GlobalErrorPercent", fix.globalErrorPercent};
        node >> Child{{"Elements", "Element"}, fix.elements};

        return node;
    }

    //! Test-facing case: the estimator Input plus the legacy golden results.
    struct LegacyCase
    {
        std::string name;
        Input input;
        double globalErrorPercent{0.0};      //!< legacy global error %
        std::vector<double> elementError{};  //!< legacy per-element relative error
        std::vector<int> refineElements{};   //!< element indices the legacy flagged to refine
    };

    inline LegacyCase loadFixture(const std::string & path)
    {
        const auto parsed = Common::loadFromFile<FixtureCase>(path, "Fixture");
        if (!parsed.has_value())
        {
            throw std::runtime_error("cannot load fixture: " + path);
        }
        const auto & fix = parsed.value();

        LegacyCase out;
        out.name = fix.name;
        out.input.targetPercent = fix.errorLimitPercent;
        out.globalErrorPercent = fix.globalErrorPercent;

        std::map<std::pair<double, double>, int> nodeIndex;
        const auto nodeOf = [&](double xpos, double ypos) -> int
        {
            const auto key = std::make_pair(xpos, ypos);
            const auto found = nodeIndex.find(key);
            if (found != nodeIndex.end())
            {
                return found->second;
            }
            const int idx = static_cast<int>(out.input.nodes.size());
            out.input.nodes.push_back(Point{xpos, ypos});
            nodeIndex.emplace(key, idx);
            return idx;
        };

        int eidx = 0;
        for (const auto & ele : fix.elements)
        {
            Element elem{};
            const double inv = 1.0 / ele.conductivity;
            elem.inverseConstitutive = {inv, 0.0, 0.0, inv};
            elem.subdomain = ele.region;

            const int vtxCount = static_cast<int>(ele.vertices.size() / 2);
            std::array<int, 4> ids = {0, 0, 0, 0};
            for (int vtx = 0; vtx < vtxCount; ++vtx)
            {
                ids[vtx] = nodeOf(ele.vertices[2 * vtx], ele.vertices[2 * vtx + 1]);
            }
            // A triangle repeats its last node so vertexIds[3] == vertexIds[2]
            // (the encoding the estimator's isTriangle relies on).
            if (vtxCount == 3)
            {
                ids[3] = ids[2];
            }
            elem.vertexIds = ids;

            for (int gpt = 0; gpt < 4; ++gpt)
            {
                elem.gaussPoints[gpt] = {ele.gaussPoints[2 * gpt], ele.gaussPoints[2 * gpt + 1]};
                elem.flux[gpt] = {ele.flux[2 * gpt], ele.flux[2 * gpt + 1]};
            }

            out.input.elements.push_back(elem);
            out.elementError.push_back(ele.legacyError);
            if (ele.refine)
            {
                out.refineElements.push_back(eidx);
            }
            ++eidx;
        }
        return out;
    }
}
