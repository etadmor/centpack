#include "centpack_2d_FD2.h"

using namespace std;

void CENTPACK::corrector_2d_FD2(doublearray3d& un, const doublearray3d& unhalf, const doublearray3d& u_star, const double& lambda, const double& mu, const double& dx, const double& dy, const double& gamma, const double& alpha, const bool& odd)
{
	long J=un.getIndex1Size() - 4;
	long K=un.getIndex2Size() - 4;
	long L=un.getIndex3Size();
	long j, k, l;
	
	doublearray1d u_star_vec(L), fluxx(L), fluxy(L);
	doublearray3d f_star(J+4,K+4,L), g_star(J+4,K+4,L);
	
	//calculate fluxes
	
	for (j=1; j<J+3; j++)
	{
		for (k=1; k<K+3; k++)
		{
			for (l=0; l<L; l++)
				u_star_vec(l)=u_star(j,k,l);
			
			flux_x(u_star_vec, fluxx, gamma);
			flux_y(u_star_vec, fluxy, gamma);
		  
			for (l=0; l<L; l++)
			{
				f_star(j,k,l)=fluxx(l);
				g_star(j,k,l)=fluxy(l);
			}
		}
	}

	if(odd)
	{
		for (j=1; j<J+2; j++)
		{
			for (k=1; k<K+2; k++)
			{
				for (l=0; l<L; l++)
					 un(j,k,l) = unhalf(j,k,l) - .5*(lambda*((f_star(j+1,k,l)-f_star(j,k,l)) + (f_star(j+1,k+1,l)-f_star(j,k+1,l))) + mu*((g_star(j,k+1,l)-g_star(j,k,l)) + (g_star(j+1,k+1,l)-g_star(j+1,k,l))));
			}
		}
	}
	
	else
	{
		for (j=2; j<J+3; j++)
		{
			for (k=2; k<K+3; k++)
			{
				for (l=0; l<L; l++)
				   un(j,k,l)= unhalf(j,k,l) - .5*(lambda*((f_star(j,k-1,l)-f_star(j-1,k-1,l)) + (f_star(j,k,l)-f_star(j-1,k,l))) + mu*((g_star(j-1,k,l)-g_star(j-1,k-1,l)) + (g_star(j,k,l)-g_star(j,k-1,l))));
			}
		}
	}
}
