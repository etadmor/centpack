#include "centpack_1d_SD3.h"

using namespace std;

void CENTPACK::C_flux_1d_SD3(const doublearray2d& u_E, const doublearray2d& u_W, const double& lambda, const double& gamma, const double& B1, doublearray2d& C)
{
	long j, l, J, L;
	
	J=C.getIndex1Size() - 4;
	L=C.getIndex2Size();
	
	doublearray1d u_ejm1(L);
	doublearray1d u_e(L);
	doublearray1d u_w(L);
	doublearray1d u_wjp1(L);
	doublearray1d H_halfp(L), H_halfm(L);

	for (j=2; j<J+2; j++)
	{
		for (l=0; l<L; l++)
		{
			u_ejm1(l) = u_E(j-1,l);
			u_e(l) = u_E(j,l);
			u_w(l) = u_W(j,l);
			u_wjp1(l) = u_W(j+1,l);
		}
		
		H_flux_1d_SD3(u_w, u_ejm1, gamma, B1, H_halfm);
		H_flux_1d_SD3(u_wjp1, u_e, gamma, B1, H_halfp);
		
		for (l=0; l<L; l++)
				C(j,l)=-lambda*(H_halfp(l)-H_halfm(l));
	}
}
