#pragma once

#include <vector>
#include <functional>

#include "SquareMatrix.hxx"

namespace MoisThermFEM {

	class IBCLinear2D;

	class BoundaryConditions2D {
	public:
		explicit BoundaryConditions2D(
				const std::vector< std::reference_wrapper< const IBCLinear2D > > & t_BCs );

		BoundaryConditions2D( const BoundaryConditions2D & other );
		BoundaryConditions2D & operator=( const BoundaryConditions2D & other );

		FenestrationCommon::SquareMatrix< double > HMatrix() const;
		FenestrationCommon::Vector< double > RVector() const;

		bool isLinear() const;

	private:
		FenestrationCommon::SquareMatrix< double > m_HMatrix;
		FenestrationCommon::Vector< double > m_RVector;
		bool m_Linear;
	};

}