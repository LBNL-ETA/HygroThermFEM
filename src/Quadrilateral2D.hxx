#pragma once

#include <vector>
#include <memory>

#include "Node2D.hxx"

namespace MoisThermFEM
{
    ////////////////////////////////////////////////////////////////////////////
    ///   QuadrilateralLinearGlobal2D
    ////////////////////////////////////////////////////////////////////////////
    class QuadrilateralLinearGlobal2D
    {
    public:
        QuadrilateralLinearGlobal2D(Node2D & t_Node1,
                                    Node2D & t_Node2,
                                    Node2D & t_Node3,
                                    Node2D & t_Node4);

        QuadrilateralLinearGlobal2D(const QuadrilateralLinearGlobal2D & t_Element);

        // Coordinates of integration (Gauss) point in global coordinate system
        double xg(const size_t IntPointIndex) const;
        double yg(const size_t IntPointIndex) const;

        // Shape function derivatives in global coordinate system
        std::vector<double> DPsiDx(const size_t IntPointIndex) const;
        std::vector<double> DPsiDy(const size_t IntPointIndex) const;

        double det(const size_t IntPointIndex) const;

        std::vector<std::size_t> nodeIndexes() const;

    private:
        QuadrilateralNodes2D m_Nodes;

        ////////////////////////////////////////////////////////////////////////////
        ///   GaussPoint2DGlobal
        ////////////////////////////////////////////////////////////////////////////

        // Handles single integration point in global coordinate space (cartesian 2D)
        class GaussPoint2DGlobal
        {
        public:
            // Nodes represent global coordinates and Index is the index of integration point
            GaussPoint2DGlobal(const Node2D & t_Node1,
                               const Node2D & t_Node2,
                               const Node2D & t_Node3,
                               const Node2D & t_Node4,
                               const size_t Index);

            double xg() const;
            double yg() const;

            std::vector<double> getDPsiDx() const;
            std::vector<double> getDPsiDy() const;

            double det() const;

        private:
            size_t m_Index;   // Index of gauss point withing global element
            double m_Xg;
            double m_Yg;
            double m_JacobiDet;
            std::vector<double> m_DPsiDx;
            std::vector<double> m_DPsiDy;
        };

        std::vector<GaussPoint2DGlobal> m_GaussPoints;
    };

}   // namespace MoisThermFEM