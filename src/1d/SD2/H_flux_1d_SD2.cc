#include "centpack_1d_SD2.h"

using namespace std;
using namespace CENTPACK;

void CENTPACK::H_flux_1d_SD2(const doublearray1d& u_w, const doublearray1d& u_e, const double& gamma, const double& B1, doublearray1d& H)
{
	long L=u_e.getIndex1Size();
	long l;
	double re, rw, a;
	
	doublearray1d f_e(L), f_w(L);
	
	spectral_radius(u_e, gamma, B1, re);
	spectral_radius(u_w, gamma, B1, rw);

	if (re>rw)
		a=re;
	else
		a=rw;
	
  	flux_x(u_e, f_e, gamma, B1);
	flux_x(u_w, f_w, gamma, B1);
	
	for (l=0; l<L; l++)
		H(l)=.5*(f_w(l)+f_e(l))-.5*a*(u_w(l)-u_e(l));
}
