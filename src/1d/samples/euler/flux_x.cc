#include "centpack_1d_FD2.h"
#include "centpack_1d_SD2.h"
#include "centpack_1d_SD3.h"

using namespace std;

void CENTPACK::flux_x(const doublearray1d& u_vector, doublearray1d& y, const double& gamma, const double& B1)
{
	double rho = u_vector(0);
	double u = u_vector(1)/rho;
	double E = u_vector(2);
	double p = (gamma - 1.0)*(E - .5*rho*pow(u,2.0));
	
	y(0) = rho*u;
	y(1) = rho*pow(u,2.0) + p;
	y(2) = u*(E + p);
}
