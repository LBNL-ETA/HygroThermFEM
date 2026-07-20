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
        //! \param maxNodeIndex Maximum node index for matrix sizing.
        SquareMatrix conductanceMatrix(size_t maxNodeIndex);

        //! Creates lumped mass matrix that includes time derivative
        //!
        //! \param maxNodeIndex Maximum node index for matrix sizing.
        //! \param DTime Timestep value for which matrix will be evaluated.
        //! \return Vector of mass values in nodes.
        std::vector<double> getLumpedMass(size_t maxNodeIndex, double DTime) const;

        //! Mass matrix in full form (not lumped)
        //!
        //! \param maxNodeIndex Maximum node index for matrix sizing.
        //! \param DTime Timestep value for which matrix will be evaluated.
        //! \return Mass matrix values.
        SquareMatrix getMassMatrix(size_t maxNodeIndex, double DTime) const;

        //! Right hand side vector of the domain.
        //! \param maxNodeIndex Maximum node index for vector sizing.
        std::vector<double> RVector(size_t maxNodeIndex) const;

        //! \brief Sets the same constant volumetric source on every element.
        void setVolumetricSource(double value);

        //! \brief Sets a constant volumetric source per element, in element creation
        //! order; the vector size must match the number of elements.
        void setVolumetricSource(const std::vector<double> & perElement);

        //! \brief Assembled consistent load vector of all element volumetric sources.
        //! \param maxNodeIndex Maximum node index for vector sizing.
        std::vector<double> volumetricSourceVector(size_t maxNodeIndex) const;

        //! Vector of fluxes for the entire domain.
        //! \param maxNodeIndex Maximum node index for vector sizing.
        std::vector<NodeFlux> flux(size_t maxNodeIndex) const;

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
