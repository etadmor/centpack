#include "centpack_1d_FD2.h"
#include "centpack_1d_SD2.h"
#include "centpack_1d_SD3.h"

using namespace std;

void CENTPACK::spectral_radius(const doublearray1d& u_vector, const double& gamma, const double& B1, double& rx)
{
	double eps = 1e-6;
	double u = u_vector(0);
	double denom = u*u + B1*(1.0 - u)*(1.0 - u);
	
	
	if ((u==0) || (u==1))
		rx = eps;
	else
		rx = (2.0*B1*u*(1.0 - u))/(denom*denom);
	//rx = 2.1;
}
