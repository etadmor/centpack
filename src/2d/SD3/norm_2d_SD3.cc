////////////////////////////////////////////////////////////////////////////
//
// CentPack -- A generic numerical solver for hyperbolic conservation
//             laws and related time dependent problems
//
// Copyright (C) 2006 Jorge Balbas and Eitan Tadmor
//
// norm.cc -- a core function of CentPack
//
// Function called by evolution
//
// This funcion calculates the norm of each conserved quantity over the
// solution domain
//
// Input (passed by reference):
//
// (1) un  --  a doublearray3d type variable holding the cell 
//             averages of u
//
// (2) dx  --  a double type variable holding the space scale dx
//
// (3) dy  --  a double type variable holding the space scale dy
//
// Output (returned by reference):
//
// (1) u_norm -- a doubearray1d variable holding the l_2 of the solution 
//               variables over the solution domain
//
////////////////////////////////////////////////////////////////////////////////

#include "centpack_2d_SD3.h"

using namespace std;

void CENTPACK::norm_2d_SD3(const doublearray3d& un, doublearray1d& u_norm, const double& dx, const double& dy)
{
	long J=un.getIndex1Size() - 4;
	long K=un.getIndex2Size() - 4;
	long L=un.getIndex3Size();
	
	long j, k, l;
	
	for (l=0; l<L; l++)
	{
		u_norm(l)=0.0;
		for (j=2; j<J+2; j++)
		{
			for (k=2; k<K+2; k++)
				u_norm(l)=u_norm(l)+pow(un(j,k,l),2.0)*dx*dy;
		}
		
		u_norm(l)=sqrt(u_norm(l));
	}
}
