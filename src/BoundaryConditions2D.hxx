#pragma once

#include <memory>
#include <vector>
#include <functional>

#include "SparceSquareMatrix.hxx"
#include "FEMunique.hxx"
#include "State.hxx"
#include "IBCLine2D.hxx"

namespace MoisThermFEM {

	/// Container for all boundary conditions. Its responsibility is to create matrix and vector of
	/// all boundary conditions.
	class BoundaryConditions2D {
	public:

		BoundaryConditions2D() = default;
		BoundaryConditions2D( const BoundaryConditions2D & other ) = delete;

		FenestrationCommon::SparceSquareMatrix< double > HMatrix() const;

		/// FenestrationCommon::SquareMatrix< double >DHMatrix() const;
		std::vector< double > RVector() const;

		bool isLinear() const;

		void updateNodeValues( const std::vector< double > & values, const Property property );

	protected:
		std::vector< std::unique_ptr< IBCLinear2D > > m_BCs;
		bool m_Linear;
	};

}