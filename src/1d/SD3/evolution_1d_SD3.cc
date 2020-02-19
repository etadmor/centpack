#include "centpack_1d_SD3.h"

using namespace std;

void CENTPACK::evolution_1d_SD3(doublearray2d& un, const double& lambda, const double& dx, const double& gamma, const double& B1)
{
	long j, l, J, L;
	
	J = un.getIndex1Size() - 4;
	L = un.getIndex2Size();
	
	doublearray1d u_norm(L);
	doublearray1d ISl(J+4), ISc(J+4), ISr(J+4);
	
	doublearray2d u_E(J+4,L);
	doublearray2d u_W(J+4,L);
	doublearray2d C0(J+4,L), C1(J+4,L), C2(J+4,L);
	
	boundary_conditions(un);
	
	norm_1d_SD3(un, u_norm, dx);
	indicators_1d_SD3(un, ISl, ISc, ISr, u_norm);
	reconstruction_1d_SD3(un, u_E, u_W, ISl, ISc, ISr);
	boundary_conditions(u_E);
	boundary_conditions(u_W);
	C_flux_1d_SD3(u_E, u_W, lambda, gamma, B1, C0);
	
	for (l=0; l<L; l++)
	{
		for (j=2; j<J+2; j++)
			un(j,l) = un(j,l) + C0(j,l);
	}
	
	boundary_conditions(un);
	
	reconstruction_1d_SD3(un, u_E, u_W, ISl, ISc, ISr);
	boundary_conditions(u_E);
	boundary_conditions(u_W);
	C_flux_1d_SD3(u_E, u_W, lambda, gamma, B1, C1);
	
	for (l=0; l<L; l++)
	{
		for (j=2; j<J+2; j++)
			un(j,l)=un(j,l) + .25*(-3.0*C0(j,l) + C1(j,l));
	}
	
	boundary_conditions(un);
	
	reconstruction_1d_SD3(un, u_E, u_W, ISl, ISc, ISr);
	boundary_conditions(u_E);
	boundary_conditions(u_W);
	C_flux_1d_SD3(u_E, u_W, lambda, gamma, B1, C2);
	
	for (l=0; l<L; l++)
	{
		for (j=2; j<J+2; j++)
			un(j,l)=un(j,l) + (1.0/12.0)*(-C0(j,l) - C1(j,l) + 8.0*C2(j,l));
	}
	
	boundary_conditions(un);
}
