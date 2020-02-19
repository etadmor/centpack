////////////////////////////////////////////////////////////////////////////
//
// CentPack -- A generic numerical solver for hyperbolic conservation
//             laws and related time dependent problems
//
// Copyright (C) 2005 Jorge Balbas and Eitan Tadmor
//
// reconstruction.cc -- a core function of CentPack
//
// Function called by evolution
//
// Requires: minmod3
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
// Output (returned by reference):
//
// (1) Four doublearray3d type variables
//
//     u_N  --  holds the point values of u at the cell interfaces
//              (x_j, y_k+1/2) over the entire solution domain
//
//     u_S  --  holds the point values of u at the cell interfaces
//              (x_j, y_k-1/2) over the entire solution domain
//
//     u_E  --  holds the point values of u at the cell interfaces
//              (x_j+1/2, y_k) over the entire solution domain
//
//     u_W  --  holds the point values of u at the cell interfaces
//              (x_j-1/2, y_k) over the entire solution domain
//
////////////////////////////////////////////////////////////////////////////////

#include "centpack_2d_SD2.h"

using namespace std;

void CENTPACK::reconstruction_2d_SD2(doublearray3d& un, doublearray3d& u_N, doublearray3d& u_S, doublearray3d& u_E, doublearray3d& u_W, const double& alpha)
{

	long J=un.getIndex1Size() - 4;
	long K=un.getIndex2Size() - 4;
	long L=un.getIndex3Size();
	long j, k, l;
	
	double ux, uy;
	
	for (l=0; l<L; l++)
	{		
		for (k=2; k<K+2; k++)
		{
			for (j=2; j<J+2; j++)
			{
				ux=minmod3(alpha*(un(j+1,k,l) - un(j,k,l)), .5*(un(j+1,k,l) - un(j-1,k,l)), alpha*(un(j,k,l) - un(j-1,k,l)));
				uy=minmod3(alpha*(un(j,k+1,l) - un(j,k,l)), .5*(un(j,k+1,l) - un(j,k-1,l)), alpha*(un(j,k,l) - un(j,k-1,l)));
				
				u_N(j,k,l)=un(j,k,l) + .5*uy;
				u_S(j,k,l)=un(j,k,l) - .5*uy;
				u_E(j,k,l)=un(j,k,l) + .5*ux;
				u_W(j,k,l)=un(j,k,l) - .5*ux;
			}
		}
	}
}
