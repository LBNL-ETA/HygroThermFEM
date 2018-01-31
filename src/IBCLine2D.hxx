#pragma once

#include <vector>

#include "SquareMatrix.hxx"
#include "Node2D.hxx"

namespace MoisThermFEM {

	class LineNodes2D;
	class LineLinearLocal1D;

	// Interface class for boundary conditions.
	class IBCLinear2D {
	public:
		IBCLinear2D() = delete;
		IBCLinear2D( const Node2D & t_Node1, const Node2D & t_Node2, const bool t_Nonlinear = false );

		std::vector< std::size_t > getNodeIndexes() const;

		virtual std::vector< double > rightHandSideVector( bool Linear = true ) const = 0;
		virtual FenestrationCommon::SquareMatrix< double > matrixA( bool Linear = true ) const = 0;

		bool isNonlinear() const;

	protected:
		static std::size_t numOfIntegrationPoints();

		static double psi( const std::size_t IntegrationPointIndex, const std::size_t Index );

		LineNodes2D m_Nodes;
		double m_Determinant;
		bool m_Nonlinear;

		/// Matrix that is base for all boundary conditions. It needs to be modified for
		/// coefficients and that will depend on type of boundary conditions
		FenestrationCommon::SquareMatrix< double > m_PsiPsiMatrix;

		/// Vector that is base for all boundary conditions. It needs to be modified for
		/// coefficients and that will depend on type of boundary conditions
		std::vector< double > m_PsiVector;
		///std::vector< double > m_Rvector; // Right hand-side vector
		///FenestrationCommon::SquareMatrix< double > m_matrixA; // Left hand-side matrix

	};

}