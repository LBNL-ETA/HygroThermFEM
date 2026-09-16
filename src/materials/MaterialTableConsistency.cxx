#include "MaterialTableConsistency.hxx"

#include <optional>
#include <sstream>
#include <utility>

#include "Material.hxx"
#include "Materials.hxx"

namespace HygroThermFEM
{
    namespace
    {
        //! A table whose points all sit at one water content is a constant: it holds at
        //! every water content and has no axis to compare.
        bool spansWaterContent(const std::vector<FenestrationCommon::point> & table)
        {
            return table.size() > 1u && table.front().x != table.back().x;
        }

        std::optional<WaterContentAxisMismatch>
          mismatchOf(const IMaterial & material,
                     std::string tableName,
                     const std::vector<FenestrationCommon::point> & table)
        {
            const double sorptionMaximum{material.sorptionCurve().back().y};
            if(!spansWaterContent(table) || table.back().x == sorptionMaximum)
            {
                return std::nullopt;
            }
            return WaterContentAxisMismatch{.materialName = material.name(),
                                            .tableName = std::move(tableName),
                                            .tableMaximum = table.back().x,
                                            .sorptionMaximum = sorptionMaximum};
        }

        void appendIfPresent(std::vector<WaterContentAxisMismatch> & mismatches,
                             std::optional<WaterContentAxisMismatch> mismatch)
        {
            if(mismatch.has_value())
            {
                mismatches.push_back(std::move(mismatch.value()));
            }
        }

        //! Default stream precision: "145.987" and "146", not the six fixed decimals
        //! std::to_string would print.
        std::string formatWaterContent(const double value)
        {
            std::ostringstream stream;
            stream << value;
            return stream.str();
        }
    }   // namespace

    std::string WaterContentAxisMismatch::message() const
    {
        return "material '" + materialName + "': " + tableName + " ends at "
               + formatWaterContent(tableMaximum)
               + " kg/m3 but the moisture storage function reaches "
               + formatWaterContent(sorptionMaximum)
               + " kg/m3; beyond its last point the table is extrapolated";
    }

    std::vector<WaterContentAxisMismatch> checkWaterContentAxes(const IMaterial & material)
    {
        std::vector<WaterContentAxisMismatch> mismatches;
        if(!material.hasSorptionCurve())
        {
            return mismatches;
        }
        if(material.hasThermalConductivityMoistureAndTemperatureDependent())
        {
            const auto conductivity{material.thermalConductivityMoistureAndTemperatureDependent()};
            appendIfPresent(mismatches,
                            mismatchOf(material,
                                       "moisture-dependent thermal conductivity",
                                       conductivity.firstTable()));
        }
        if(material.hasLiquidTransportationCurve())
        {
            appendIfPresent(
              mismatches,
              mismatchOf(material, "liquid transport coefficient", material.liquidTransportationCurve()));
        }
        if(material.hasDiffusionResistanceFactorMoistureDependent())
        {
            appendIfPresent(mismatches,
                            mismatchOf(material,
                                       "moisture-dependent vapor resistance factor",
                                       material.diffusionResistanceFactorMoistureDependent()));
        }
        return mismatches;
    }

    std::vector<WaterContentAxisMismatch> checkWaterContentAxes(const Materials & materials)
    {
        std::vector<WaterContentAxisMismatch> mismatches;
        for(const auto & materialName : materials.getSolidMaterials())
        {
            auto found{checkWaterContentAxes(materials.material(materialName))};
            mismatches.insert(mismatches.end(),
                              std::make_move_iterator(found.begin()),
                              std::make_move_iterator(found.end()));
        }
        return mismatches;
    }
}   // namespace HygroThermFEM
