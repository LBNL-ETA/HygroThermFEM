#include "FEMMatrices.hxx"

#include "VectorOperators.hxx"

namespace HygroThermFEM
{
    SquareMatrix steadyStateLeftHandSide(SingleDomain & domain)
    {
        auto K = domain.m_Elements.conductanceMatrix();
        K += domain.m_BCs.HMatrix(0);

        return K;
    }

    std::vector<double> steadyStateRightHandSide(SingleDomain & domain)
    {
        return domain.m_BCs.RVector();
    }

    SquareMatrix transientM_K_H_Matrix(SingleDomain & domain, double t_DTime, size_t timestepIndex)
    {
        const auto M = domain.m_Elements.getLumpedMass(t_DTime);
        auto M_K_H = domain.m_Elements.conductanceMatrix();
        M_K_H = M_K_H.addDiagonal(M);
        M_K_H += domain.m_BCs.HMatrix(timestepIndex);

        return M_K_H;
    }

    std::vector<double> transientMT_R_Vector(SingleDomain & domain,
                                             const std::vector<double> & t_PreviousTimestepValues,
                                             double t_DTime,
                                             size_t timestepIndex)
    {
        const auto M{domain.m_Elements.getLumpedMass(t_DTime)};
        const auto R{domain.m_BCs.RVector(timestepIndex) + domain.m_Elements.RVector()};

        return t_PreviousTimestepValues * M + R;
    }
}   // namespace HygroThermFEM