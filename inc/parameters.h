#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <cmath>

namespace parameter
{

	const double dt = 1e-7;
	const double a = 1;
	const double A_0 = 0.6;
	const double K_a = 0.5;
	
	const double l_min = 0.005*std::sqrt(A_0);
	const double l_new = 0.01*std::sqrt(A_0);
	const double A_min = 0.5*A_0;
	const double A_max = 2.0*A_0;
	
	const double Lambda = 0.1;
	const double Gamma = 0.5;

};

#endif // PARAMETERS_H