#pragma once

#include <memory>

#include "Node2D.hxx"
#include "Quadrilateral2D.hxx"
#include "SparceSquareMatrix.hxx"
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
		IQLEMatrix2D( const QuadrilateralLinearGlobal2D & t_Element );

		// Integrate matrix over all points of integration
		virtual FenestrationCommon::SparceSquareMatrix< double >
		integrate( const std::vector< double > & t_Values ) const final;

	protected:
		virtual void calculateMatrixInIntegrationPoint(
				const std::vector< double > & t_Values,
				const std::size_t t_IntegrationPointIndex,
				FenestrationCommon::SparceSquareMatrix< double > & t_Matrix ) const final;

		const QuadrilateralLinearGlobal2D & m_Global2D;

		std::vector< FenestrationCommon::SparceSquareMatrix< double > > m_IntegrationMatrix;

	};

	//////////////////////////////////////////////////////////////////////////////
	///  QLEConductance2D
	//////////////////////////////////////////////////////////////////////////////

	// Class to handle conductance matrix in global coordinate system
	class QLEConductance2D : public IQLEMatrix2D {
	public:
		QLEConductance2D( const QuadrilateralLinearGlobal2D & t_Element );

	};

	//////////////////////////////////////////////////////////////////////////////
	///  QLEConductanceDerivative2D
	//////////////////////////////////////////////////////////////////////////////

	// Handles conductance part with derivative term
	class QLEConductanceDerivative2D : public IQLEMatrix2D {
	public:
		QLEConductanceDerivative2D( const QuadrilateralLinearGlobal2D & t_Element );

		// This updates integration matrix with new derivative values
		void updateIntegrationMatrix( const std::vector< double > & t_Values );

		void clearIntegrationMatrix();

	};

	//////////////////////////////////////////////////////////////////////////////
	///  QLECapacitance2D
	//////////////////////////////////////////////////////////////////////////////

	// Class to handle capacitance matrix in global coordinate system
	class QLECapacitance2D : public IQLEMatrix2D {
	public:
		QLECapacitance2D( const QuadrilateralLinearGlobal2D & t_Element );

	};

	/// Keeping function pointers for QLEConductanceDerivative2D in Elements array
	struct DerivativeFunction {
		DerivativeFunction( const std::shared_ptr< IValue > & fixedTerm,
												const std::shared_ptr< IValue > & derivativeTerm );

		std::shared_ptr< MoisThermFEM::IValue > fixedTerm;
		std::shared_ptr< MoisThermFEM::IValue > derivativeTerm;
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
											const Node2D & t_Node4, const Material & t_Material );

		FenestrationCommon::SparceSquareMatrix< double > conductanceMatrix() const;

		FenestrationCommon::SparceSquareMatrix< double > conductanceDerivativeMatrix();

		FenestrationCommon::SparceSquareMatrix< double > capacitanceMatrix() const;

		Node2D & getNode( std::size_t index );

		bool haveBothNodes( const Node2D & t_Node1, const Node2D & t_Node2 ) const;

		std::vector< std::size_t > nodeIndexes() const;

		const Material & getMaterial() const;

	protected:
		/// TODO: This did not work with reference_wrapper and it should. Check later.
		/// Reminder: Introduce pair of curve pointer and Property so that curve knows what to use
		std::vector< std::shared_ptr< IValue > > m_Conductance;
		std::vector< std::shared_ptr< IValue > > m_Capacitance;
		std::vector< DerivativeFunction > m_DerivativeConductance;

		const Material & m_Material;

	private:
		std::vector< Node2D > m_Node;

		QuadrilateralLinearGlobal2D m_Global2D;
		QLECapacitance2D m_QLECapacitance2D;
		QLEConductance2D m_QLEConductance2D;
		/// This one depends on functions and must be stored for every DerivativeConductance submatrix
		std::vector< QLEConductanceDerivative2D > m_QLEDerivativeConductance;
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
