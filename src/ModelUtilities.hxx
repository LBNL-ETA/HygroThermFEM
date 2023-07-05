#pragma once

#include <string>
#include <memory>
#include <functional>

#include "../src/GravityVector.hxx"

namespace HygroThermFEM
{
    struct MultiDomain;
    struct SingleDomain;
    class IElementLinear2D;

    using ElementFactory =
      std::function<std::unique_ptr<IElementLinear2D>(const size_t index1,
                                                      const size_t index2,
                                                      const size_t index3,
                                                      const size_t index4,
                                                      const std::string & materialName)>;

    //! Adds element into domain
    //!< Node 1 index
    //!< Node 2 index
    //!< Node 3 index
    //!< Node 4 index
    //!< SolidMaterial that will be assigned to the element
    void createElement(SingleDomain & domain,
                       size_t index1,
                       size_t index2,
                       size_t index3,
                       size_t index4,
                       const std::string & materialName);

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

    //! \brief Checks if domain equations are linear
    bool isLinear(SingleDomain & domain);

    //! @brief Deletes Geometry, boundary conditions and materials
    void clearModel(SingleDomain & domain);

    //! @brief Deletes Geometry, boundary conditions and materials
    void clearModel(MultiDomain & domain);

    //! \brief Sets new gravity vector and performs new calculations
    //!
    //! @param gravityVector Direction of gravity
    void setGravityVector(MultiDomain & domain,
                          const FenestrationCommon::GravityVector & gravityVector);
}   // namespace HygroThermFEM