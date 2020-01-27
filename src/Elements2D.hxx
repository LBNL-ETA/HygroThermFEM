#pragma once

#include <memory>
#include "Element2D.hxx"
#include "SquareMatrix.hxx"

namespace HygroThermFEM
{
    //! \brief Container class to hold all elements connected into global matrix. This is only for
    //! elements and not for boundary conditions
    class ElementsLinear2D
    {
    public:
        explicit ElementsLinear2D() = default;

        //! Conductance matrix calculated for the entire domain.
        SquareMatrix conductanceMatrix();

        //! Creates lumped mass matrix that includes time derivative
        //!
        //! \param DTime Timestep value for which matrix will be evaluated.
        //! \return Vector of mass values in nodes.
        std::vector<double> getLumpedMass(double DTime);

        //! Mass matrix in full form (not lumped)
        //!
        //! \param DTime Timestep value for which matrix will be evaluated.
        //! \return Mass matrix values.
        SquareMatrix getMassMatrix(double DTime);

        //! Right hand side vector of the domain.
        std::vector<double> RVector() const;

        //! Vector of fluxes for the entire domain.
        std::vector<NodeFlux> flux() const;

        //! Flag that shows if domain can be solved as linear.
        bool isLinear() const;

        //! \brief Returns the element with given nodes. Note that one two elements can contain same
        //! edge, however, that edge will have have nodes in different order and becuase of that
        //! node order passed to this function is important.
        //!
        //! \param index1 First node of the segment.
        //! \param index2 Second node of the segment.
        //! \return Returns element (if exists)
        IElementLinear2D * findElement(size_t index1, size_t index2);

        //! \brief Assigns new element to the domain.
        //!
        //! \param el New element that will be moved into domain.
        void assignElement(std::unique_ptr<IElementLinear2D> && el);

        //! \brief Returns all elements in the domain.
        const std::vector<std::unique_ptr<IElementLinear2D>> & elements() const;

        //! \brief Clears elements to be ready for next simulation
        void clearElements();

    protected:
        std::vector<std::unique_ptr<IElementLinear2D>> m_Elements;
    };

}   // namespace HygroThermFEM
