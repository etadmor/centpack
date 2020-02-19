#include "centpack_1d_FD2.h"
#include "centpack_1d_SD2.h"
#include "centpack_1d_SD3.h"

using namespace std;

void CENTPACK::spectral_radius(const doublearray1d& u, const double& gamma, const double& B1, double& rx)
{
	using namespace std;

	double rho, vx, p, A, B, cf;
	
	rho = u(0);
	vx = u(1)/rho;
	p = u(6)-.5*((pow(u(1),2.0) + pow(u(2),2.0) + pow(u(3),2.0))/rho) - .5* (pow(u(4),2.0) + pow(u(5),2.0) + pow(B1,2.0));
	A = 2.0*p/rho;
	B = (pow(B1,2.0) + pow(u(4),2.0) + pow(u(5),2.0))/rho;
	cf = sqrt(.5*(A+B+sqrt(pow(A+B,2.0) - 4.0*A*pow(B1,2.0)/rho)));
	rx = fabs(vx)+cf;
}
