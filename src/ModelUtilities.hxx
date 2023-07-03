#pragma once

#include <string>
#include "../src/GravityVector.hxx"

namespace HygroThermFEM
{
    struct MultiDomain;

    //! \brief Creates element with material reference
    //! @param index1 Node 1 index
    //! @param index2 Node 2 index
    //! @param index3 Node 3 index
    //! @param index4 Node 4 index
    //! @param materialName Material name assigned to the element
    void createElement(MultiDomain & domain,
                       size_t index1,
                       size_t index2,
                       size_t index3,
                       size_t index4,
                       const std::string & materialName);

    //! \brief Sets new gravity vector and performs new calculations
    //!
    //! @param gravityVector Direction of gravity
    void setGravityVector(MultiDomain & domain,
                          const FenestrationCommon::GravityVector & gravityVector);

    //! @brief Deletes Geometry and boundary conditions
    void clearModel(MultiDomain & domain);
}   // namespace HygroThermFEM