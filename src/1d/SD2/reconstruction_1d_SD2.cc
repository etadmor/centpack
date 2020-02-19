#include "centpack_1d_SD2.h"

using namespace std;
using namespace CENTPACK;

void CENTPACK::reconstruction_1d_SD2(doublearray2d& un, doublearray2d& u_E, doublearray2d& u_W, const double& alpha)
{
	long J=un.getIndex1Size() - 4;
	long L=un.getIndex2Size();
	long j, l;

	double u, s;

	for (l=0; l<L; l++)
	{
		for (j=2; j<J+2; j++)
		{
			u=un(j,l);
			s = minmod3(alpha*(un(j+1,l)-un(j,l)), .5*(un(j+1,l)-un(j-1,l)), alpha*(un(j,l)-un(j-1,l)));
		
			u_E(j,l) = u + .5*s;
			u_W(j,l) = u - .5*s;
		}
	}
}
