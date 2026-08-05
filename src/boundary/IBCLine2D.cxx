#include <cmath>
#include <stdexcept>

#include "IBCLine2D.hxx"
#include "IntegrationPoints.hxx"
#include "LineLocal1D.hxx"
#include "Elements2D.hxx"
#include "Nodes.hxx"

namespace HygroThermFEM
{
    IBCLinear2D::IBCLinear2D(Nodes & nodePool,
                             const size_t index1,
                             const size_t index2,
                             const bool t_Linear) :
        m_Nodes(nodePool.getNode(index1), nodePool.getNode(index2)),
        m_Linear(t_Linear)
    {
        m_Determinant =
          0.5
          * sqrt(pow(m_Nodes[0].X() - m_Nodes[1].X(), 2) + pow(m_Nodes[0].Y() - m_Nodes[1].Y(), 2));

        // Create psi matrix and vector
        for(std::size_t i = 0; i < numOfIntegrationPoints(); ++i)
        {
            for(std::size_t j = 0; j < numOfBCNodes; ++j)
            {
                for(std::size_t k = 0; k < numOfBCNodes; ++k)
                {
                    m_PsiPsiMatrix(j, k) += m_Determinant * psi(i, j) * psi(i, k);
                }
                m_PsiVector[j] += m_Determinant * psi(i, j);
            }
        }
    }

    double IBCLinear2D::gaussPointProperty(const std::size_t integrationPointIndex,
                                           const Variable variable) const
    {
        double result{0.0};
        for(std::size_t idx = 0; idx < numOfBCNodes; ++idx)
        {
            result += psi(integrationPointIndex, idx) * m_Nodes[idx].property(variable);
        }
        return result;
    }

    std::array<std::size_t, numOfBCNodes> IBCLinear2D::getNodeIndexes() const
    {
        return {m_Nodes[0].getNodeNumber(), m_Nodes[1].getNodeNumber()};
    }

    size_t IBCLinear2D::numOfIntegrationPoints()
    {
        return IntegrationPoints2D::Instance().count1D();
    }

    double IBCLinear2D::psi(const size_t IntegrationPointIndex, const size_t Index)
    {
        return LineLinearLocal1D::Instance().Psi(IntegrationPointIndex, Index);
    }

    bool IBCLinear2D::isLinear() const
    {
        return m_Linear;
    }

    INode2D & IBCLinear2D::getNode(const std::size_t index) const
    {
        if(index > numOfBCNodes)
        {
            throw std::runtime_error("Index out of range.");
        }
        return m_Nodes[index];
    }

}   // namespace HygroThermFEM
