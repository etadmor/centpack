#include "centpack_1d_FD2.h"
#include "centpack_1d_SD2.h"
#include "centpack_1d_SD3.h"

using namespace std;

void CENTPACK::flux_x(const doublearray1d& u_vector, doublearray1d& y, const double& gamma, const double& B1)
{
	double u = u_vector(0);
	
	y(0) = .5*u*u;
}
