#include "centpack_1d_SD2.h"

using namespace std;
using namespace CENTPACK;

void CENTPACK::evolution_1d_SD2(doublearray2d& un, const double& lambda, const double& dx, const double& gamma, const double& alpha, const double& B1)
{
	long j, l, J, L;
	
	J = un.getIndex1Size() - 4;
	L = un.getIndex2Size();
	
	doublearray2d u_E(J+4,L);
	doublearray2d u_W(J+4,L);
	doublearray2d C0(J+4,L), C1(J+4,L);
	
	boundary_conditions(un);
	
	reconstruction_1d_SD2(un, u_E, u_W, alpha);
	boundary_conditions(u_E);
	boundary_conditions(u_W);
	C_flux_1d_SD2(u_E, u_W, lambda, gamma, B1, C0);
	
	for (l=0; l<L; l++)
	{
		for (j=2; j<J+2; j++)
			un(j,l) = un(j,l) + C0(j,l);
	}
	
	boundary_conditions(un);
	
	reconstruction_1d_SD2(un, u_E, u_W, alpha);
	boundary_conditions(u_E);
	boundary_conditions(u_W);
	C_flux_1d_SD2(u_E, u_W, lambda, gamma, B1, C1);
	
	for (l=0; l<L; l++)
	{
		for (j=2; j<J+2; j++)
			un(j,l) = un(j,l) + .5*(C1(j,l)-C0(j,l));
	}
	
	boundary_conditions(un);
}
