#include "centpack_1d_FD2.h"

using namespace std;

void CENTPACK::corrector_1d_FD2(doublearray2d& un, const doublearray2d& uj_half, const doublearray2d& un_half, const double& lambda, const double& gamma, const double& B1, const bool& odd)
{
	long j, l, J, L;
	
	J=un.getIndex1Size() - 4;
	L=un.getIndex2Size();
	
	doublearray1d u_vector(L);
	doublearray1d f_vector(L);
	doublearray2d f_half(J+4,L);
	
	for (j=1; j<J+3; j++)
	{
		for (l=0; l<L; l++)
			u_vector(l) = un_half(j,l);
		
		flux_x(u_vector, f_vector, gamma, B1);
		
		for (l=0; l<L; l++)
			f_half(j,l) = f_vector(l);
	}
	
	if (odd)
	{
		for (j=1; j<J+2; j++)
		{
			for (l=0; l<L; l++)
				un(j,l) = uj_half(j,l) - lambda*(f_half(j+1,l) - f_half(j,l));
		}
	}
	
	else
	{
		for (j=2; j<J+3; j++)
		{
			for (l=0; l<L; l++)
				un(j,l) = uj_half(j,l) - lambda*(f_half(j,l) - f_half(j-1,l));
		}
	}
}
