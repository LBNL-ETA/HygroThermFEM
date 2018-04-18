#pragma once

#include <vector>

#include "SquareMatrix.hxx"
#include "Node2D.hxx"
#include "Vector.hxx"
#include "Material.hxx"

namespace MoisThermFEM {

	const std::size_t numOfBCNodes = 2;

	class LineNodes2D;
	class LineLinearLocal1D;

	// Interface class for boundary conditions in linear 2D finite elementsCreator.
	class IBCLinear2D {
	public:
		virtual ~IBCLinear2D();

		IBCLinear2D() = delete;
		IBCLinear2D( const Node2D & t_Node1, const Node2D & t_Node2, const bool t_Linear = true );

		std::vector< std::size_t > getNodeIndexes() const;
		Node2D & getNode( const std::size_t index );

		/// Every boundary condition will return matrices and vectors. There are two type of matrices.
		/// First one is called H matrix and will include sum of all matrices that do not have first
		/// derivative of coefficients. Second one is delta H matrix that is based on coefficients
		/// derivative. Reason for keeping them separate is that matrix equation for non-linear case
		/// need to know what matrices are coming from coefficients derivative, because they will be on
		/// left hand side, next to unknown.
		virtual FenestrationCommon::Vector< double > R_Vector() const = 0;
		virtual FenestrationCommon::SquareMatrix< double > H_Matrix() const = 0;

		/// DHMatrix seems unnecessary because solution will converge anyway.
		/// virtual FenestrationCommon::SquareMatrix< double > D_HMatrix() const;

		bool isLinear() const;

	protected:

		static std::size_t numOfIntegrationPoints();

		static double psi( const std::size_t IntegrationPointIndex, const std::size_t Index );

		double getIntegratedProperty( const Property t_Property ) const;
		double getIntegratedDeltaProperty( const Property t_Property ) const;

		LineNodes2D m_Nodes;
		double m_Determinant;
		bool m_Linear;

		/// Matrix that is base for all boundary conditions. It needs to be modified for
		/// coefficients and that will depend on type of boundary conditions
		FenestrationCommon::SquareMatrix< double > m_PsiPsiMatrix;

		/// Vector that is base for all boundary conditions. It needs to be modified for
		/// coefficients and that will depend on type of boundary conditions
		FenestrationCommon::Vector< double > m_PsiVector;

	};

}