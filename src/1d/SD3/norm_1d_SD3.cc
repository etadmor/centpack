#include "centpack_1d_SD3.h"

using namespace std;

void CENTPACK::norm_1d_SD3(const doublearray2d& un, doublearray1d& u_norm, const double& dx)
{
	long J=un.getIndex1Size() - 4;
	long L=un.getIndex2Size();
	
	long j, l;
	
	for (l=0; l<L; l++)
	{
		u_norm(l)=0.0;
		for (j=2; j<J+2; j++)
			u_norm(l)=u_norm(l) + pow(un(j,l),2.0)*dx;
		
		u_norm(l) = sqrt(u_norm(l));
	}
}
