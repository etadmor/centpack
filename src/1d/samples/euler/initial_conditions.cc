#include "centpack_1d_FD2.h"
#include "centpack_1d_SD2.h"
#include "centpack_1d_SD3.h"

using namespace std;

void CENTPACK::initial_conditions(doublearray2d& un, const double& x_init, const double& x_final, const double& gamma, doublearray1d& x)
{
	
	long Nu = un.getIndex1Size() - 4;
	long L = un.getIndex2Size();
	long J = L;
	long nu, j, l;
	long nu1 = Nu/2 + 2;
	double dx;
	double p_left;
	double p_right;
	
	doublearray1d u_left(L), u_right(L);
	
	dx = (x_final - x_init)/Nu;
	
	x(0) = x_init - 1.5*dx;
	
	for (nu = 1; nu < Nu+4; nu++)
		x(nu) = x(nu-1) + dx;
	
	p_left = 1.0;
	p_right = 0.1;
	
	u_left(0) = 1.0;
	u_left(1) = 0.0;
	u_left(2) = p_left/(gamma - 1.0);
	
	u_right(0) = 0.125; 
	u_right(1) = 0.0;
	u_right(2) = p_right/(gamma - 1.0);
	
	for(l = 0; l < L; l++)
	{
		for (nu = 0; nu < nu1; nu++)
			un(nu,l) = u_left(l);
		
		for (nu = nu1; nu < Nu + 4; nu++)
			un(nu,l) = u_right(l);
	}
}
