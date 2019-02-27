#pragma once

#include <memory>
#include "Element2D.hxx"
#include "SquareMatrix.hxx"

namespace HygroThermFEM
{
    /// Container class to hold all elements connected into global matrix. This is only for elements
    /// and not for boundary conditions
    class ElementsLinear2D
    {
    public:
        explicit ElementsLinear2D() = default;

        FenestrationCommon::SquareMatrix conductanceMatrix();

        /// Creates lumped mass matrix that includes time derivative
        std::vector<double> getLumpedMass(double DTime);
        FenestrationCommon::SquareMatrix getMassMatrix(double DTime);

        std::vector<double> RVector() const;

        std::vector<NodeFlux> flux() const;

        bool isLinear() const;

        void updateNodeValues(const std::vector<double> & values,
                              BaseVariable property,
                              bool updatePreviousValue = true);

        IElementLinear2D * findElement(size_t index1, size_t index2);

        void assignElement(std::unique_ptr<IElementLinear2D> && el);

        const std::vector<std::unique_ptr<IElementLinear2D>> & elements() const;

    protected:
        std::vector<std::unique_ptr<IElementLinear2D>> m_Elements;
    };

}   // namespace HygroThermFEM
