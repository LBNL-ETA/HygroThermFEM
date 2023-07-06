#include <memory>
#include <stdexcept>

#include "LineLocal1D.hxx"

#include "IntegrationPoints.hxx"
#include "Node2D.hxx"

namespace HygroThermFEM
{
    ////////////////////////////////////////////////////////////////////////////
    //   LineElement1DLinearLocal::LineLinearLocalShapeFunctions1D
    ////////////////////////////////////////////////////////////////////////////
    LineLinearLocal1D::LineLinearLocalShapeFunctions1D::LineLinearLocalShapeFunctions1D(
      LocalPoint1D const & t_Point) :
        ILocalShapeFunctions1DLine(t_Point)
    {
        m_Psi.push_back(0.5 * (1 - t_Point.ksi));
        m_Psi.push_back(0.5 * (1 + t_Point.ksi));
    }

    ////////////////////////////////////////////////////////////////////////////
    //   LineElement1DLinearLocal
    ////////////////////////////////////////////////////////////////////////////

    LineLinearLocal1D & LineLinearLocal1D::Instance()
    {
        static LineLinearLocal1D m_Instance;
        return m_Instance;
    }

    LineLinearLocal1D::LineLinearLocal1D()
    {
        std::vector<LocalPoint1D> aPoints = IntegrationPoints2D::Instance().getPoints1D();
        for(const auto & point : aPoints)
        {
            m_Ksi.push_back(std::unique_ptr<LineLinearLocalShapeFunctions1D>(
              std::make_unique<LineLinearLocalShapeFunctions1D>(
                LineLinearLocalShapeFunctions1D(point))));
        }
    }

    double LineLinearLocal1D::Psi(const size_t IntegrationPointIndex, const size_t Index)
    {
        if(IntegrationPointIndex >= m_Ksi.size())
        {
            throw std::runtime_error(
              "Integration point index out of range. Routine LineElement1DLinearLocal::Psi.");
        }

        return m_Ksi[IntegrationPointIndex]->Psi(Index);
    }

}   // namespace HygroThermFEM
