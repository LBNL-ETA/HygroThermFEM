#pragma once

#include <cmath>

namespace MoisThermFEM {
	
	const double ConvergenceError = 1e-5;
	const size_t MaxIterations = 1000;

}

namespace Constants {
	static const double PI = atan(1)*4;
	static const double STEFANBOLTZMANN = 5.6697E-8;

	//// TODO: Keep these constant for now. Gases can calculate these properties.
	static const double Density_Air = 1.2922;
	static const double Cp_Air = 1;
	static const double K_Air = 0.025;

	//// TODO: Check if water and ice properties need to be recalculated
	static const double Density_Water = 1000;
	static const double Cp_Water = 4184;
	static const double K_Water = 0.591;

	static const double Density_Ice = 916.7;
	static const double Cp_Ice = 2108;
	static const double K_Ice = 2.22;

}
