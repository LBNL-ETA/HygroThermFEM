#pragma once

#include "Node2D.hxx"
#include "Quadrilateral2D.hxx"
#include "SquareMatrix.hxx"

namespace MoisThermFEM {	

	// Constant that holds number of nodes in certain elements
	const std::size_t numOfQadrilateralNodes = 4;

	//////////////////////////////////////////////////////////////////////////////
	//  IQuadrilateralElement2D
	//////////////////////////////////////////////////////////////////////////////

	class IElementQuadrilateral2D {
	public:
		virtual ~IElementQuadrilateral2D() = default;

		IElementQuadrilateral2D( const Node2D & t_Node1, const Node2D & t_Node2, const Node2D & t_Node3,
		                 const Node2D & t_Node4 );

		std::vector< size_t > nodeIndexes() const;

	protected:
		QuadrilateralLinearGlobal2D m_Element;
		QuadrilateralNodes2D m_ElementNodes;
	};

	//////////////////////////////////////////////////////////////////////////////
	//  IQLEMatrix2D
	//////////////////////////////////////////////////////////////////////////////

	// Abstract class that forces users to perform matrix calculation in inherited version
	// Depending on equations, matrices will have different calculation methods (for example
	// capacitance and conductance matrices have different form)
	class IQLEMatrix2D {
	public:
		IQLEMatrix2D( const double t_Value );

		FenestrationCommon::SquareMatrix< double > getMatrix() const;

		void calculate();

	protected:
		virtual FenestrationCommon::SquareMatrix< double > calculateMatrixInIntegrationPoint( 
			const std::size_t t_IntegrationPointIndex ) const = 0;

		FenestrationCommon::SquareMatrix< double > m_Matrix;
		const double m_Value;
	};

	//////////////////////////////////////////////////////////////////////////////
	//  IQLEConductance2D
	//////////////////////////////////////////////////////////////////////////////

	// Class to handle conductance matrix in global coordinate system
	class IQLEConductance2D : public IElementQuadrilateral2D, public IQLEMatrix2D {
	public:
		IQLEConductance2D( const Node2D & t_Node1, const Node2D & t_Node2, const Node2D & t_Node3,
		                 const Node2D & t_Node4, const double t_Value );

	protected:
		FenestrationCommon::SquareMatrix< double > calculateMatrixInIntegrationPoint(
			const std::size_t t_IntegrationPointIndex ) const override;

	};

	//////////////////////////////////////////////////////////////////////////////
	//  IQLECapacitance2D
	//////////////////////////////////////////////////////////////////////////////

	class IQLECapacitance2D : public IElementQuadrilateral2D, public IQLEMatrix2D {
	public:
		IQLECapacitance2D( const Node2D & t_Node1, const Node2D & t_Node2, const Node2D & t_Node3,
		                 const Node2D & t_Node4, const double t_Value );

	protected:
		FenestrationCommon::SquareMatrix< double > calculateMatrixInIntegrationPoint(
			const size_t t_IntegrationPointIndex ) const override;

	};

	//////////////////////////////////////////////////////////////////////////////
	//  ElementThermalLinear2D
	//////////////////////////////////////////////////////////////////////////////

	// Handles linear 2D element (4 nodes)
	class ElementThermalLinear2D {
	public:
		ElementThermalLinear2D( const Node2D & t_Node1, const Node2D & t_Node2, const Node2D & t_Node3,
		                 const Node2D & t_Node4, const double t_Cond, const double t_Rho = 1,
		                 const double t_Cp = 1 );

		ElementThermalLinear2D( const ElementThermalLinear2D & t_Element );

		std::vector< size_t > nodeIndexes() const;

		FenestrationCommon::SquareMatrix< double > conductivity() const;
		FenestrationCommon::SquareMatrix< double > rhoCp() const;

	private:
		IQLECapacitance2D m_Capacitance;
		IQLEConductance2D m_Conductance;
	};

}
