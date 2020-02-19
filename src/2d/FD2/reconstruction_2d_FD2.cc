//////////////////////////////////////////////////////////////////////////////////
//
// reconstruction.cc -- a core function of
//
// Central Hyperbolic Solver (C) Jorge Balbas and Eitan Tadmor, Nov. 2004
//
// Function called by evolution
//
// This function implements a 3rd order Kurganov-Levy CWENO reconstruction
// of the interface values of the solution u from its cell averages along
// the x- and z-directions
// 
// Input (passed by reference):
//
// (1) un  --  a doublearray3d type variable holding the cell averages of u
//             at time t=t^n over the discretized solution domain
//
// (2) Six doublearray2d arrays holding the value of the smoothness 
//     indicators used to calculate the non-linear weights of the 
//     non-oscillatory reconstruction along the x- and z-directions:
//
// Output (returned by reference):
//
// (1) Four doublearray3d type variables
//
//
////////////////////////////////////////////////////////////////////////////////

#include "centpack_2d_FD2.h"

using namespace std;

void CENTPACK::reconstruction_2d_FD2(const doublearray3d& un, doublearray3d& unhalf, const double& alpha, const bool& odd)
{

	long J=un.getIndex1Size() - 4;
	long K=un.getIndex2Size() - 4;
	long L=un.getIndex3Size();
	long j, k, l;
	
	doublearray3d u_prime(J+4,K+4,L);
	doublearray3d u_qrime(J+4,K+4,L);
	
	for (l=0; l<L; l++)
	{		
		for (k=1; k<K+3; k++)
		{
			for (j=1; j<J+3; j++)
			{
				u_prime(j,k,l)=minmod3(alpha*(un(j+1,k,l) - un(j,k,l)), .5*(un(j+1,k,l) - un(j-1,k,l)), alpha*(un(j,k,l) - un(j-1,k,l)));
				u_qrime(j,k,l)=minmod3(alpha*(un(j,k+1,l) - un(j,k,l)), .5*(un(j,k+1,l) - un(j,k-1,l)), alpha*(un(j,k,l) - un(j,k-1,l)));
			}
		}
		
		if (odd)
		{
			for (k=1; k<K+2; k++)
			{
				for (j=1; j<J+2; j++)
					unhalf(j,k,l) = .25*((un(j,k,l) + un(j+1,k,l) + un(j,k+1,l) + un(j+1,k+1,l)) + .25*((u_prime(j,k,l) - u_prime(j+1,k,l)) + (u_prime(j,k+1,l) - u_prime(j+1,k+1,l)) + (u_qrime(j,k,l) - u_qrime(j,k+1,l)) + (u_qrime(j+1,k,l) - u_qrime(j+1,k+1,l))));
			}
		}
		
		else
		{
			for (k=2; k<K+3; k++)
			{
				for (j=2; j<J+3; j++)
					unhalf(j,k,l)= .25*((un(j,k-1,l) + un(j-1,k-1,l) + un(j,k,l) + un(j-1,k,l)) + .25*((u_prime(j-1,k-1,l) - u_prime(j,k-1,l)) + (u_prime(j-1,k,l) - u_prime(j,k,l)) + (u_qrime(j-1,k-1,l) - u_qrime(j-1,k,l)) + (u_qrime(j,k-1,l)-u_qrime(j,k,l))));
			}
		}
	}
	
	//boundary_conditions(unhalf);
}
