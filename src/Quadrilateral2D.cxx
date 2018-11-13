#include <cmath>
#include <stdexcept>
#include <cassert>

#include "Quadrilateral2D.hxx"
#include "IntegrationPoints.hxx"
#include "QuadrilateralLocal2D.hxx"

namespace MoisThermFEM
{
    ////////////////////////////////////////////////////////////////////////////
    //   QuadrilateralLinearGlobal2D
    ////////////////////////////////////////////////////////////////////////////

    QuadrilateralLinearGlobal2D::QuadrilateralLinearGlobal2D(Node2D & t_Node1,
                                                             Node2D & t_Node2,
                                                             Node2D & t_Node3,
                                                             Node2D & t_Node4) :
        m_Nodes(t_Node1, t_Node2, t_Node3, t_Node4)
    {
        // Setting up gauss points in global coordinate system
        auto & aElement = IntegrationPoints2D::Instance();
        auto numOfPoints = aElement.count2D();
        for(unsigned i = 0; i < numOfPoints; ++i)
        {
            m_GaussPoints.emplace_back(t_Node1, t_Node2, t_Node3, t_Node4, i);
        }
    }

    QuadrilateralLinearGlobal2D::QuadrilateralLinearGlobal2D(
      const QuadrilateralLinearGlobal2D & t_Element) :
        m_Nodes(t_Element.m_Nodes)
    {
        m_GaussPoints.clear();
        for(auto const & aPoint : t_Element.m_GaussPoints)
        {
            m_GaussPoints.push_back(aPoint);
        }
    }

    double QuadrilateralLinearGlobal2D::xg(const size_t IntPointIndex) const
    {
        assert(IntPointIndex < m_GaussPoints.size());
        return m_GaussPoints[IntPointIndex].xg();
    }

    double QuadrilateralLinearGlobal2D::yg(const size_t IntPointIndex) const
    {
        assert(IntPointIndex < m_GaussPoints.size());
        return m_GaussPoints[IntPointIndex].yg();
    }

    std::vector<double> QuadrilateralLinearGlobal2D::DPsiDx(const size_t IntPointIndex) const
    {
        assert(IntPointIndex < m_GaussPoints.size());
        return m_GaussPoints[IntPointIndex].getDPsiDx();
    }

    std::vector<double> QuadrilateralLinearGlobal2D::DPsiDy(const size_t IntPointIndex) const
    {
        assert(IntPointIndex < m_GaussPoints.size());
        return m_GaussPoints[IntPointIndex].getDPsiDy();
    }

    double QuadrilateralLinearGlobal2D::det(const size_t IntPointIndex) const
    {
        assert(IntPointIndex < m_GaussPoints.size());
        return m_GaussPoints[IntPointIndex].det();
    }

    std::vector<std::size_t> QuadrilateralLinearGlobal2D::nodeIndexes() const
    {
        return m_Nodes.getNodeIndexes();
    }

    ////////////////////////////////////////////////////////////////////////////
    //   GaussPoint2DGlobal
    ////////////////////////////////////////////////////////////////////////////

    QuadrilateralLinearGlobal2D::GaussPoint2DGlobal::GaussPoint2DGlobal(const Node2D & t_Node1,
                                                                        const Node2D & t_Node2,
                                                                        const Node2D & t_Node3,
                                                                        const Node2D & t_Node4,
                                                                        const size_t Index) :
        m_Index(Index),
        m_Xg(0),
        m_Yg(0)
    {
        std::vector<std::reference_wrapper<const Node2D>> node = {
          t_Node1, t_Node2, t_Node3, t_Node4};

        auto & aElement = QuadrilateralLinearLocal2D::Instance();

        // Calculate jacobian matrix manually rather then using ublas from boost. Matrix is
        // fixed to 2x2 and this should be the fastest way
        double j11 = 0;
        double j12 = 0;
        double j21 = 0;
        double j22 = 0;

        // global derivatives
        double dxdksi = 0;
        double dydksi = 0;
        double dxdeta = 0;
        double dydeta = 0;

        for(unsigned i = 0; i < 4; ++i)
        {
            auto aNode = node[i].get();
            m_Xg += aNode.X() * aElement.Psi(Index, i);
            m_Yg += aNode.Y() * aElement.Psi(Index, i);

            dxdksi += aNode.X() * aElement.PsiDKsi(Index, i);
            dydksi += aNode.Y() * aElement.PsiDKsi(Index, i);
            dxdeta += aNode.X() * aElement.PsiDEta(Index, i);
            dydeta += aNode.Y() * aElement.PsiDEta(Index, i);

            j11 += aNode.X() * aElement.PsiDKsi(Index, i);
            j12 += aNode.Y() * aElement.PsiDKsi(Index, i);
            j21 += aNode.X() * aElement.PsiDEta(Index, i);
            j22 += aNode.Y() * aElement.PsiDEta(Index, i);
        }

        m_JacobiDet = j11 * j22 - j21 * j12;

        // integrate inverse matrix now
        auto temp = j11;
        j11 = j22 / m_JacobiDet;
        j22 = temp / m_JacobiDet;
        j12 = -j12 / m_JacobiDet;
        j21 = -j21 / m_JacobiDet;

        for(size_t i = 0; i < 4; ++i)
        {
            m_DPsiDx.push_back(j11 * aElement.PsiDKsi(Index, i) + j12 * aElement.PsiDEta(Index, i));
            m_DPsiDy.push_back(j21 * aElement.PsiDKsi(Index, i) + j22 * aElement.PsiDEta(Index, i));
        }
    }

    double QuadrilateralLinearGlobal2D::GaussPoint2DGlobal::xg() const
    {
        return m_Xg;
    }

    double QuadrilateralLinearGlobal2D::GaussPoint2DGlobal::yg() const
    {
        return m_Yg;
    }

    std::vector<double> QuadrilateralLinearGlobal2D::GaussPoint2DGlobal::getDPsiDx() const
    {
        return m_DPsiDx;
    }

    std::vector<double> QuadrilateralLinearGlobal2D::GaussPoint2DGlobal::getDPsiDy() const
    {
        return m_DPsiDy;
    }

    double QuadrilateralLinearGlobal2D::GaussPoint2DGlobal::det() const
    {
        return m_JacobiDet;
    }

}   // namespace MoisThermFEM
