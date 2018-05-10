#pragma once

#include "Elements2D.hxx"
#include "BoundaryCondition2DFactory.hxx"
#include "Element2DFactory.hxx"

namespace MoisThermFEM {

	/// Class to keep elements and boundary conditions together. One domain per Thermal, Mass and
	/// Pressure governing equations.
	class Domain {
	public:
		explicit Domain( const Property property );

		/// Calculates steady state solution
		std::vector< double > steadyState();

		/// Calculates next timestep value from current values
		std::vector< double > transient( std::vector< double > & currentStateValues, double t_DTime );

		void createConvectionBC( const Node2D & t_Node1, const Node2D & t_Node2,
								double t_ConvectionCoefficient, double t_AirTemperature );

		void createTemperatureBC( Node2D & t_Node1, Node2D & t_Node2, double t_Temp1, double t_Temp2 );

		void createTemperatureBC( Node2D & t_Node1, Node2D & t_Node2, const double t_Temp );

		void createFluxBC( Node2D & t_Node1, Node2D & t_Node2, const double t_Flux );

		void createBlackBodyRadiationBC( const Node2D & t_Node1, const Node2D & t_Node2,
																		 const double t_Emissivity,
																		 const double t_RadiationTemperature );

		void createMoistureBC( const Node2D & t_Node1, const Node2D & t_Node2,
															 const double t_ConvectiveCoefficient, const double t_AirHumidity,
															 const double t_AirTemperature );

		void createThermalElement( const Node2D & t_Node1, const Node2D & t_Node2, const Node2D & t_Node3,
															 const Node2D & t_Node4, const Material & mat );

		void createMoistureElement( const Node2D & t_Node1, const Node2D & t_Node2, const Node2D & t_Node3,
																const Node2D & t_Node4, const Material & mat );

		IElementLinear2D* findElement( const Node2D & t_Node1, const Node2D & t_Node2 );

	protected:

		friend class MultiDomain;

		void updateNodeValues( const std::vector< double > & values, const Property property );

		FenestrationCommon::SparceSquareMatrix< double > steadyStateLeftHandSide();
		std::vector< double > steadyStateRightHandSide() const;

		/// In matrix equations some structures are showing up in both (linear and nonlinear) cases
		/// and those matrix operations are separated into functions.
		/// This function retrieves M+K+H matrix
		FenestrationCommon::SparceSquareMatrix< double > transientM_K_H_Matrix( const double t_DTime );
		/// FenestrationCommon::SquareMatrix< double > transientDH_Matrix();

		/// This function retrieves M*U+R vector (where U is state variable)
		std::vector< double >
		transientMT_R_Vector( std::vector< double > & t_PreviousSolution, const double t_DTime );

		bool isLinear() const;

		Element2DFactory m_Elements;
		BoundaryCondition2DFactory m_BCs;
		Property m_Property;
	};

}