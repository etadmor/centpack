#include "centpack_2d_FD2.h"

using namespace std;

void CENTPACK::predictor_2d_FD2(const doublearray3d& un, doublearray3d& u_star, const double& lambda, const double& mu, const double& dx, const double& dy, const double& gamma, const double& alpha)
{
	long J=un.getIndex1Size() - 4;
	long K=un.getIndex2Size() - 4;
	long L=un.getIndex3Size();
	long j, k, l;

	double fu, fv, fw;
	double gu, gv, gw;
	
	doublearray1d u_vector(L), fluxx(L), fluxy(L);
	doublearray3d f(J+4,K+4,L), g(J+4,K+4,L);
	doublearray3d f_prime(J+4,K+4,L), g_qrime(J+4,K+4,L);
	
	//calculate fluxes
	
	for (j=0; j<J+4; j++)
	{
		for (k=0; k<K+4; k++)
		{
			for (l=0; l<L; l++)
				u_vector(l)=un(j,k,l);
			
			flux_x(u_vector, fluxx, gamma);
			flux_y(u_vector, fluxy, gamma);
		  
			for (l=0; l<L; l++)
			{
				f(j,k,l)=fluxx(l);
				g(j,k,l)=fluxy(l);
			}
		}
	}

//calculate flux derivatives

  	for (l=0; l<L; l++)
	{
		for (k=1; k<K+3; k++)
		{
			for (j=1; j<J+3; j++)
			{
				fu=alpha*(f(j+1,k,l)-f(j,k,l));
				fv=.5*(f(j+1,k,l)-f(j-1,k,l));
				fw=alpha*(f(j,k,l)-f(j-1,k,l));
				f_prime(j,k,l)=minmod3(fu, fv, fw);
				
				gu = alpha*(g(j,k+1,l)-g(j,k,l));
				gv = .5*(g(j,k+1,l)-g(j,k-1,l));
				gw = alpha*(g(j,k,l)-g(j,k-1,l));
				g_qrime(j,k,l)=minmod3(gu, gv, gw);
			}
		}
	}

//calculate predicted values

  	for(l=0; l<L; l++)
	{
		for (j=1; j<J+3; j++)
		{
			for (k=1; k<K+3; k++)
				u_star(j,k,l)=un(j,k,l) - .5*(lambda*f_prime(j,k,l) + mu*g_qrime(j,k,l));
		}
	}
}
