#pragma once

#include <vector>
#include <memory>

#include "SquareMatrix.hxx"

namespace MoisThermFEM {

	class ILineLinear2D;

	class BoundaryConditions2D {
	public:
		explicit BoundaryConditions2D( const std::vector< std::unique_ptr< ILineLinear2D > > & t_BCs );

		BoundaryConditions2D( const BoundaryConditions2D & other );

		BoundaryConditions2D & operator=( const BoundaryConditions2D & other );

		FenestrationCommon::SquareMatrix< double > matrixA() const;

		std::vector< double > vectorR() const;

	private:
		FenestrationCommon::SquareMatrix< double > m_MatrixA;
		std::vector< double > m_vectorR;
	};

}