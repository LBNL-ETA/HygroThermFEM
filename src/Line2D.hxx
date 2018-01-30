#pragma once

#include <vector>

#include "SquareMatrix.hxx"
#include "Node2D.hxx"

namespace MoisThermFEM {

	class LineNodes2D;
	class LineLinearLocal1D;

	// Interface class for boundary conditions.
	class ILineLinear2D {
	public:
		ILineLinear2D() = delete;
		ILineLinear2D( const Node2D & t_Node1, const Node2D & t_Node2 );

		std::vector< std::size_t > getNodeIndexes() const;

		std::vector< double > rightHandSideVector() const;
		FenestrationCommon::SquareMatrix< double > matrixA() const;

	protected:
		static std::size_t numOfIntegrationPoints();

		static double psi( const std::size_t IntegrationPointIndex, const std::size_t Index );

		LineNodes2D m_Nodes;
		double m_Determinant;
		std::vector< double > m_Rvector; // Right hand-side vector
		FenestrationCommon::SquareMatrix< double > m_matrixA; // Left hand-side matrix

	};

}