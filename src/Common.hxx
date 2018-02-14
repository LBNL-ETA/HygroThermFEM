#pragma once

namespace MoisThermFEM {
	
	const double ConvergenceError = 1e-5;
	const size_t MaxIterations = 1000;

}

namespace Constants {
	static const double STEFANBOLTZMANN = 5.6697E-8;


	//// TODO: Keep these constant for now. Gases can calculate these properties.
	static const double Density_AIR = 1.2922;
	static const double Cp_Air = 1;

	//// TODO: Check if water properties need to be recalculated
	static const double Water_Density = 1000;
	static const double Ice_Density = 916.7;
}
