#include "centpack_1d_FD2.h"

using namespace std;

void CENTPACK::reconstruction_1d_FD2(doublearray2d& un, doublearray2d& uj_half, const double& alpha, const bool& odd)
{

	using namespace std;

	long J=un.getIndex1Size() - 4;
	long L=un.getIndex2Size();
	long j, l;
	
	doublearray2d u_prime(J+4,L);
	
	for (j=1; j<J+3; j++)
	{
		for (l=0; l<L; l++)
		{
			u_prime(j,l) = minmod3(alpha*(un(j,l) - un(j-1,l)), .5*(un(j+1,l) - un(j-1,l)), alpha*(un(j+1,l) - un(j,l)));
		}
	}
	
	if (odd)
	{
		for (j=1; j<J+2; j++)
		{
			for (l=0; l<L; l++)
				uj_half(j,l) = .5*(un(j+1,l) + un(j,l)) + .125*(u_prime(j,l) - u_prime(j+1,l));
		}
	}
	
	else
	{
		for (j=2; j<J+3; j++)
		{
			for (l=0; l<L; l++)
				uj_half(j,l) = .5*(un(j,l) + un(j-1,l)) + .125*(u_prime(j-1,l) - u_prime(j,l));
		}
	}
}
