#include <stdexcept>

#include "LocalShapeFunctions.hxx"
#include "Node2D.hxx"

namespace MoisThermFEM
{
    //////////////////////////////////////////////////////////////////////////////////////////////
    //  ILocalShapeFunctions2DQuadrilateral
    //////////////////////////////////////////////////////////////////////////////////////////////

    ILocalShapeFunctions2DQuadrilateral::ILocalShapeFunctions2DQuadrilateral(
      const LocalPoint2D & t_Point) :
        m_Ksi(t_Point.ksi),
        m_Eta(t_Point.eta)
    {}

    double ILocalShapeFunctions2DQuadrilateral::Psi(const size_t Index) const
    {
        if(Index >= m_Psi.size())
        {
            throw std::runtime_error(
              "Index out of range. Routine ILocalShapeFunctions2DQuadrilateral::Psi.");
        }
        return m_Psi[Index];
    }

    double ILocalShapeFunctions2DQuadrilateral::PsiDKsi(const size_t Index) const
    {
        if(Index >= m_PsiDKsi.size())
        {
            throw std::runtime_error(
              "Index out of range. Routine ILocalShapeFunctions2DQuadrilateral::PsiDKsi.");
        }
        return m_PsiDKsi[Index];
    }

    double ILocalShapeFunctions2DQuadrilateral::PsiDEta(const size_t Index) const
    {
        if(Index >= m_PsiDEta.size())
        {
            throw std::runtime_error(
              "Index out of range. Routine ILocalShapeFunctions2DQuadrilateral::PsiDEta.");
        }
        return m_PsiDEta[Index];
    }

    const std::vector<double> & ILocalShapeFunctions2DQuadrilateral::VPsi() const
    {
        return m_Psi;
    }

    const std::vector<double> & ILocalShapeFunctions2DQuadrilateral::VPsiDKsi() const
    {
        return m_PsiDKsi;
    }

    const std::vector<double> & ILocalShapeFunctions2DQuadrilateral::VPsiDEta() const
    {
        return m_PsiDEta;
    }

    //////////////////////////////////////////////////////////////////////////////////////////////
    //  ILocalShapeFunctions1DLine
    //////////////////////////////////////////////////////////////////////////////////////////////

    ILocalShapeFunctions1DLine::ILocalShapeFunctions1DLine(const LocalPoint1D & t_Point) :
        m_Ksi(t_Point.ksi)
    {}

    double ILocalShapeFunctions1DLine::Psi(const size_t Index) const
    {
        if(Index >= m_Psi.size())
        {
            throw std::runtime_error(
              "Index out of range. Routine ILocalShapeFunctions1DLine::Psi.");
        }
        return m_Psi[Index];
    }

}   // namespace MoisThermFEM
