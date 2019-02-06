#pragma once

#include <vector>
#include <memory>

#include "Node2D.hxx"

namespace HygroThermFEM
{
    ////////////////////////////////////////////////////////////////////////////
    ///   QuadrilateralLinearGlobal2D
    ////////////////////////////////////////////////////////////////////////////

    //! \brief Calculates element determinant, shape function derivatives, global coordinates
    //! of element's Gauss points
    //!
    //! These calculations are necessary in every differential function integration.
    class QuadrilateralLinearGlobal2D
    {
    public:
        //! Constructor of the element
        QuadrilateralLinearGlobal2D(Node2D & t_Node1,   //!< Element's node 1
                                    Node2D & t_Node2,   //!< Element's node 2
                                    Node2D & t_Node3,   //!< Element's node 3
                                    Node2D & t_Node4    //!< Element's node 4
        );

        QuadrilateralLinearGlobal2D(const QuadrilateralLinearGlobal2D & t_Element);

        //! x coordinate of Gauss integration point
        double xg(const size_t IntPointIndex   //! Gauss integration point index
                  ) const;

        //! y coordinate of Gauss integration point
        double yg(const size_t IntPointIndex   //! Gauss integration point index
                  ) const;

        //! Shape function derivative over x in global coordinate system
		const std::vector< double > & DPsiDx(
			const size_t IntPointIndex   //! Gauss integration point index
		) const;

        //! Shape function derivative over y in global coordinate system
		const std::vector< double > & DPsiDy(
			const size_t IntPointIndex   //! Gauss integration point index
		) const;

        //! Element determinant in Gauss point of integration
        double det(const size_t IntPointIndex   //! Gauss integration point index
                   ) const;

        std::vector<std::size_t> nodeIndexes() const;

    private:
        //! Simple storage for nodes. This storage is adopted to always have four nodes.
        QuadrilateralNodes2D m_Nodes;

        ////////////////////////////////////////////////////////////////////////////
        ///   GaussPoint2DGlobal
        ////////////////////////////////////////////////////////////////////////////

        //! \brief Helper class that is used to calculate shape functions and it derivatives in
        //! global coordinate system for single integration point
        class GaussPoint2DGlobal
        {
        public:
            // Nodes represent global coordinates and Index is the index of integration point
            GaussPoint2DGlobal(const Node2D & t_Node1, //!< Node 1 of the element
                               const Node2D & t_Node2, //!< Node 2 of the element
                               const Node2D & t_Node3, //!< Node 3 of the element
                               const Node2D & t_Node4, //!< Node 4 of the element
                               const size_t Index //!< Integration point index
                               );

            //! x coordinate of integration point
            double xg() const;

            //! y coordinate of integration point
            double yg() const;

            //! Shape function derivatives over x coordinate
            const std::vector<double> & getDPsiDx() const;

            //! Shape function derivatives over y coordinate
            const std::vector<double> & getDPsiDy() const;

            //! Element determinant
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

}   // namespace HygroThermFEM
