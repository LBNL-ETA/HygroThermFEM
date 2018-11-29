#pragma once

#include <vector>

#include "SquareMatrix.hxx"
#include "Node2D.hxx"
#include "Material.hxx"

namespace MoisThermFEM
{
    const std::size_t numOfBCNodes = 2;

    class LineNodes2D;
    class LineLinearLocal1D;

    //! \brief Interface class for boundary conditions in linear 2D finite elementsCreator.
    //!
    //! Any new boundary condition created in MoisThermFEM engine must be inherited from this class.
    //! It provides necessary interface for rest of the program to be able to communicate with
    //! boundary condition. In finite element world, every boundary condition will return matrix
    //! and vector. Matrix is what stand next to state variable (unknown) that program is seeking
    //! solution.
    class IBCLinear2D
    {
    public:
        virtual ~IBCLinear2D() = default;
        IBCLinear2D(const IBCLinear2D & ibc) = default;

        IBCLinear2D() = delete;

        //! Basic construction of boundary condition.
        IBCLinear2D(size_t index1,         //!< Index of first node.
                    size_t index2,         //!< Index of second node.
                    bool t_Linear = true   //!< Linearity of the boundary condition.
                    //!< This variable should be set to false if boundary condition introduce
                    //!< non-linearity into domain.
        );

        //! Returns node indexes
        std::vector<std::size_t> getNodeIndexes() const;

        //! Returns node with given index.
        Node2D & getNode(std::size_t index) const;

        //! Virtual function that will calculate vector from the boundary condition.
        //! Any inherited boundary condition must define its own vector calculations.
        virtual std::vector<double> R_Vector() const = 0;

        //! Virtual function that will calculate matrix from the boundary condition.
        //! Any inherited boundary condition must define its own matrix calculations.
        virtual FenestrationCommon::SquareMatrix H_Matrix() const = 0;

        //! Shows if boundary condition is linear or not.
        bool isLinear() const;

    protected:
        //! Returns number of integration points.
        static std::size_t numOfIntegrationPoints();

        //! Calculates psi value for given integration point and shape function
        static double psi(std::size_t IntegrationPointIndex,   //!< Integration point index.
                          std::size_t Index                    //!< Shape function index
        );

        LineNodes2D m_Nodes;
        bool m_Linear;

        //! Matrix that is base for all boundary conditions. It needs to be modified for
        //! coefficients and that will depend on type of boundary conditions
        FenestrationCommon::SquareMatrix m_PsiPsiMatrix;

        //! Vector that is base for all boundary conditions (it depends on psi functions).
        //! It needs to be modified for coefficients and that will depend on type of boundary
        //! conditions
        std::vector<double> m_PsiVector;
    };

}   // namespace MoisThermFEM