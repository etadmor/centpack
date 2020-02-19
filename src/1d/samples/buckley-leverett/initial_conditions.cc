#include "centpack_1d_FD2.h"
#include "centpack_1d_SD2.h"
#include "centpack_1d_SD3.h"

using namespace std;

void CENTPACK::initial_conditions(doublearray2d& un, const double& x_init, const double& x_final, const double& gamma, doublearray1d& x)
{
	
	long J=un.getIndex1Size() - 4;
	long L=un.getIndex2Size();
	long j, l;
	long j1 = J/2 + 2;
	double dx;
	
	dx=(x_final - x_init)/J;
	
	x(0) = x_init - 2.0*dx;
	
	for (j=1; j<J+4; j++)
		x(j) = x(j-1) + dx;
	
	for(l=0; l<L; l++)
	{
		for (j=0; j<j1; j++)
			un(j,l) = 1.0;
		for (j=j1; j<J+4; j++)
			un(j,l) = 0.0;
	}
}
