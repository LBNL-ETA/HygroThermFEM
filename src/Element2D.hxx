#pragma once

#include <memory>

#include "Node2D.hxx"
#include "Quadrilateral2D.hxx"
#include "SquareMatrix.hxx"
#include "Material.hxx"
#include "Functions.hxx"

namespace MoisThermFEM {

	// Constant that holds number of nodes in certain elementsCreator
	const std::size_t numOfQuadrilateralNodes = 4;
	// const std::size_t numOfIntegrationPoints = 4;

	//////////////////////////////////////////////////////////////////////////////
	///  IQLEMatrix2D
	//////////////////////////////////////////////////////////////////////////////

	// Abstract class that forces users to perform matrix calculation in inherited version
	// Depending on equations, matrices will have different calculation methods (for example
	// capacitance and conductance matrices have different form)
	class IQLEMatrix2D {
	public:
		IQLEMatrix2D( const QuadrilateralLinearGlobal2D & t_Element, const double t_Value1,
									const double t_Value2,
									const double t_Value3, const double t_Value4 );

		FenestrationCommon::SquareMatrix< double > getMatrix() const;

		// Integrate matrix over all points of integration
		void integrate();

	protected:
		virtual FenestrationCommon::SquareMatrix< double > calculateMatrixInIntegrationPoint(
				const std::size_t t_IntegrationPointIndex ) const = 0;

		FenestrationCommon::SquareMatrix< double > m_Matrix;
		const std::vector< double > m_Values;

		const QuadrilateralLinearGlobal2D & m_Global2D;

		std::vector< FenestrationCommon::SquareMatrix< double > > m_IntegrationMatrix;

	};

	//////////////////////////////////////////////////////////////////////////////
	///  QLEConductance2D
	//////////////////////////////////////////////////////////////////////////////

	// Class to handle conductance matrix in global coordinate system
	class QLEConductance2D : public IQLEMatrix2D {
	public:
		QLEConductance2D(
				const QuadrilateralLinearGlobal2D & t_Element, const double t_Value1, const double t_Value2,
				const double t_Value3,
				const double t_Value4 );

	protected:
		FenestrationCommon::SquareMatrix< double > calculateMatrixInIntegrationPoint(
				const std::size_t t_IntegrationPointIndex ) const override;

	};

	//////////////////////////////////////////////////////////////////////////////
	///  QLECapacitance2D
	//////////////////////////////////////////////////////////////////////////////

	// Class to handle capacitance matrix in global coordinate system
	class QLECapacitance2D : public IQLEMatrix2D {
	public:
		QLECapacitance2D( const QuadrilateralLinearGlobal2D & t_Element, const double t_Value1,
											const double t_Value2,
											const double t_Value3,
											const double t_Value4 );

	protected:
		FenestrationCommon::SquareMatrix< double > calculateMatrixInIntegrationPoint(
				const std::size_t t_IntegrationPointIndex ) const override;

	};

	//////////////////////////////////////////////////////////////////////////////
	///  IElementLinear2D
	//////////////////////////////////////////////////////////////////////////////

	/// Class that handles creation of conductance and capacitance matrices in linear
	/// 2D world. This class will be used by multiple governing equations since
	/// basis of matrix creation are identical with only difference in what coefficients
	/// are passed
	class IElementLinear2D {
	public:
		IElementLinear2D( const Node2D & t_Node1, const Node2D & t_Node2, const Node2D & t_Node3,
											const Node2D & t_Node4 );

		FenestrationCommon::SquareMatrix< double > conductanceMatrix() const;

		FenestrationCommon::SquareMatrix< double > capacitanceMatrix() const;

		Node2D & getNode( const std::size_t index );

		std::vector< std::size_t > nodeIndexes() const;

	protected:
		std::vector< Node2D > m_Node;

		/// TODO: This did not work with reference_wrapper and it should. Check later.
		/// Reminder: Introduce pair of curve pointer and Property so that curve knows what to use
		std::vector< std::shared_ptr< MoisThermFEM::IValue > > m_Conductance;
		std::vector< std::shared_ptr< MoisThermFEM::IValue > > m_Capacitance;

	private:
		QuadrilateralLinearGlobal2D m_Global2D;
	};

	//////////////////////////////////////////////////////////////////////////////
	///  ElementThermalLinear2D
	//////////////////////////////////////////////////////////////////////////////

	class ElementThermalLinear2D : public IElementLinear2D {
	public:
		ElementThermalLinear2D( const Node2D & t_Node1, const Node2D & t_Node2, const Node2D & t_Node3,
														const Node2D & t_Node4, const Material & mat );

	};

	//////////////////////////////////////////////////////////////////////////////
	///  ElementMoistureLinear2D
	//////////////////////////////////////////////////////////////////////////////

	class ElementMoistureLinear2D : public IElementLinear2D {
	public:
		ElementMoistureLinear2D( const Node2D & t_Node1, const Node2D & t_Node2, const Node2D & t_Node3,
														 const Node2D & t_Node4, const Material & mat );

	};

}
