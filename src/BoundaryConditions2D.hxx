#pragma once

#include <vector>
#include <functional>

#include "SquareMatrix.hxx"

namespace MoisThermFEM {

	class IBCLinear2D;

	/// Container for all boundary conditions. Its responsibility is to create matrix and vector of
	/// all boundary conditions.
	class BoundaryConditions2D {
	public:
		explicit BoundaryConditions2D(
				std::vector< std::reference_wrapper< IBCLinear2D > > & t_BCs );

		FenestrationCommon::SquareMatrix< double > HMatrix() const;
		FenestrationCommon::SquareMatrix< double >DHMatrix() const;
		FenestrationCommon::Vector< double > RVector() const;

		bool isLinear() const;

		void updateNodeTemperatures( const std::vector< double > & temperatures );

	private:
		std::vector< std::reference_wrapper< IBCLinear2D > > m_BCs;
		bool m_Linear;
	};

}