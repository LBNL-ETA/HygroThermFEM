#pragma once

#include <vector>

#include "MaterialMissingProperties.hxx"
#include "Common.hxx"

namespace HygroThermFEM
{
    struct MultiDomain;

    //! \brief Checks validity of materials for any simulation
    //!
    //! \param simulationType Enumerator for simulation type (SteadyState or Transient)
    //! \return Information about materials with all missing properties.
    [[nodiscard]] MaterialsErrorCheckVector
      checkForMaterialsValidity(const MultiDomain & domain, SimulationType simulationType);
}   // namespace HygroThermFEM