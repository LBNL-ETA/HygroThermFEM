#pragma once

#include <memory>
#include <vector>
#include <functional>

#include "SquareMatrix.hxx"
#include "FEMunique.hxx"
#include "State.hxx"
#include "IBCLine2D.hxx"

namespace MoisThermFEM
{
    //! \brief Container for all boundary conditions.
    //!
    //! Beside simple container that will hold all kinds of boundary conditions, responsibility of
    //! this class is to create assembled vector and matrix for all boundary conditions.
    class BoundaryConditions2D
    {
    public:
        //! Simple constructor
        explicit BoundaryConditions2D();

        //! Making copy of boundary conditions container is not allowed. Since all changes
        //! and updates should be kept at one place.
        BoundaryConditions2D(const BoundaryConditions2D & other) = delete;

        //! Assembled H matrix from all boundary conditions.
        FenestrationCommon::SquareMatrix HMatrix() const;

        //! Assembled R vector from all boundary conditions.
        std::vector<double> RVector() const;

        //! Returns linearity of the problem. It is important for domain to know what set of
        //! equations should be applied.
        bool isLinear() const;

        //! Used to update new values to nodes once solver finds new solution for any state
        //! variable.
        void updateNodeValues(const std::vector<double> & values, const BaseVariable property,
                              bool updatePreviousValue = true);

        //! Assign new boundary condition to the pool
        void assignBC(std::unique_ptr<IBCLinear2D> && bc);

    protected:
        std::vector<std::unique_ptr<IBCLinear2D>> m_BCs;
        bool m_Linear;
    };

}   // namespace MoisThermFEM