#pragma once

#include "Node2D.hxx"
#include "Quadrilateral2D.hxx"
#include "SquareMatrix.hxx"

namespace MoisThermFEM {	

	// Constant that holds number of nodes in certain elements
	const std::size_t numOfQuadrilateralNodes = 4;

	//////////////////////////////////////////////////////////////////////////////
	///  IQuadrilateralElement2D
	//////////////////////////////////////////////////////////////////////////////

	class IElementQuadrilateral2D {
	public:
		virtual ~IElementQuadrilateral2D() = default;

		IElementQuadrilateral2D( const Node2D & t_Node1, const Node2D & t_Node2, const Node2D & t_Node3,
		                 const Node2D & t_Node4 );

		std::vector< std::size_t > nodeIndexes() const;

	protected:
		QuadrilateralLinearGlobal2D m_Element;
		QuadrilateralNodes2D m_ElementNodes;
	};

	//////////////////////////////////////////////////////////////////////////////
	///  IQLEMatrix2D
	//////////////////////////////////////////////////////////////////////////////

	// Abstract class that forces users to perform matrix calculation in inherited version
	// Depending on equations, matrices will have different calculation methods (for example
	// capacitance and conductance matrices have different form)
	class IQLEMatrix2D {
	public:
		IQLEMatrix2D( const double t_Value );

		FenestrationCommon::SquareMatrix< double > getMatrix() const;

		// Integrate matrix over all points of integration
		void integrate();

	protected:
		virtual FenestrationCommon::SquareMatrix< double > calculateMatrixInIntegrationPoint( 
			const std::size_t t_IntegrationPointIndex ) const = 0;

		FenestrationCommon::SquareMatrix< double > m_Matrix;
		const double m_Value;
	};

	//////////////////////////////////////////////////////////////////////////////
	///  QLEConductance2D
	//////////////////////////////////////////////////////////////////////////////

	// Class to handle conductance matrix in global coordinate system
	class QLEConductance2D : public IElementQuadrilateral2D, public IQLEMatrix2D {
	public:
		QLEConductance2D( const Node2D & t_Node1, const Node2D & t_Node2, const Node2D & t_Node3,
		                 const Node2D & t_Node4, const double t_Value );

	protected:
		FenestrationCommon::SquareMatrix< double > calculateMatrixInIntegrationPoint(
			const std::size_t t_IntegrationPointIndex ) const override;

	};

	//////////////////////////////////////////////////////////////////////////////
	///  QLECapacitance2D
	//////////////////////////////////////////////////////////////////////////////

	// Class to handle capacitance matrix in global coordinate system
	class QLECapacitance2D : public IElementQuadrilateral2D, public IQLEMatrix2D {
	public:
		QLECapacitance2D( const Node2D & t_Node1, const Node2D & t_Node2, const Node2D & t_Node3,
		                 const Node2D & t_Node4, const double t_Value );

	protected:
		FenestrationCommon::SquareMatrix< double > calculateMatrixInIntegrationPoint(
			const std::size_t t_IntegrationPointIndex ) const override;

	};

	//////////////////////////////////////////////////////////////////////////////
	///  IElementlLinear2D
	//////////////////////////////////////////////////////////////////////////////

	/// Class that handles creation of conductance and capacitance matrices in linear
	/// 2D world. This class will be inhertied by multple governing equations since
	/// basis of matrix creation are identical with only difference in what coefficients
	/// are passed
	class IElementLinear2D {
	public:
		IElementLinear2D( const Node2D & t_Node1, const Node2D & t_Node2, const Node2D & t_Node3,
											const Node2D & t_Node4, const double t_Conductance, const double t_Capacitance );

		std::vector< std::size_t > nodeIndexes() const;

		FenestrationCommon::SquareMatrix< double > conductanceMatrix() const;
		FenestrationCommon::SquareMatrix< double > capacitanceMatrix() const;

	private:
		QLEConductance2D m_Conductance;
		QLECapacitance2D m_Capacitance;
	};

	//////////////////////////////////////////////////////////////////////////////
	///  ElementThermalLinear2D
	//////////////////////////////////////////////////////////////////////////////

	// Handles linear 2D element (4 nodes)
	class ElementThermalLinear2D : public IElementLinear2D {
	public:
		ElementThermalLinear2D( const Node2D & t_Node1, const Node2D & t_Node2, const Node2D & t_Node3,
		                 const Node2D & t_Node4, const double t_Cond = 1, const double t_Rho = 1,
		                 const double t_Cp = 1 );

	};

}
