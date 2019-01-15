#pragma once

#include <memory>
#include <functional>
#include "Element2D.hxx"
#include "SquareMatrix.hxx"

namespace MoisThermFEM
{
    /// Container class to hold all elements connected into global matrix. This is only for elements
    /// and not for boundary conditions
    class ElementsLinear2D
    {
    public:
        explicit ElementsLinear2D() = default;

        FenestrationCommon::SquareMatrix conductanceMatrix();

        /// Creates lumped mass matrix that includes time derivative
        std::vector<double> getLumpedMass(const double DTime);
        FenestrationCommon::SquareMatrix getMassMatrix(const double DTime);

        std::vector<double> RVector() const;

        std::vector<NodeFlux> flux() const;

        bool isLinear() const;

        void updateNodeValues(const std::vector<double> & values, const BaseVariable property,
                              bool updatePreviousValue = true);

        IElementLinear2D * findElement(const size_t index1, const size_t index2);

        void assignElement(std::unique_ptr<IElementLinear2D> && el);

    protected:
        /// FenestrationCommon::SquareMatrix< double > m_Conductance;
        /// FenestrationCommon::SquareMatrix< double > m_Capacitance;

        std::vector<std::unique_ptr<IElementLinear2D>> m_Elements;
    };

}   // namespace MoisThermFEM