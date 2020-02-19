////////////////////////////////////////////////////////////////////////////
//
// CentPack -- A generic numerical solver for hyperbolic conservation
//             laws and related time dependent problems
//
// Copyright (C) 2006 Jorge Balbas and Eitan Tadmor
//
// indicators.cc -- a core function of CentPack
//
// Function called by evolution
//
// Requires: sum (see prototype below)
//
// This function calculates the smoothness indicators used to calculate the 
// non-linear weights of the non-oscillatory reconstruction of u as the 
// normed-average of the smoothness indicators of its components
// 
// Input (passed by reference):
// 
// (1) un  --  a doublearray3d type variable holding the cell averages of u
//              at time t=t^n over the discretized solution domain
//
// 
// (3) u_norm  --  doablearray1d holding the values of the l_2 norm of cell 
//                 averages of each conserved quantity over the solution domain
//
// Output (returned by reference):
//
// (1) Six doublearray2d arrays holding the value of the smoothness 
//     indicators along the x- and z-directions:
//     
//     ISl  --  left reconstruction
//
//     IScx  --  central reconstruction along x-direction
//
//     ISr  --  right reconstruction
//
//     ISb  --  bottom reconstruction
//
//     IScy  --  central reconstruction along y-direction
//
//     ISt  --  top reconstruction
//
////////////////////////////////////////////////////////////////////////////////

#include "centpack_2d_SD3.h"

using namespace std;

void CENTPACK::indicators_2d_SD3(doublearray3d& un, doublearray2d& ISl, doublearray2d& IScx, doublearray2d& ISr, doublearray2d& ISb, doublearray2d& IScy, doublearray2d& ISt, const doublearray1d& u_norm)
{
	
	long J=un.getIndex1Size() - 4;
	long K=un.getIndex2Size() - 4;
	long L=un.getIndex3Size();
	long j, k, l;
	
	doublearray1d pl(2), pcx(3), pr(2);
	doublearray1d pt(2), pcy(3), pb(2);
	doublearray1d ISl_vector(L), IScx_vector(L), ISr_vector(L), ISb_vector(L), IScy_vector(L), ISt_vector(L);
	
	double eps=0.000001;
	
	for (k=2; k<K+2; k++)
	{
		for (j=2; j<J+2; j++)
		{
			for (l=0; l<L; l++)
			{
				pl(0)=un(j,k,l);
				pl(1)=un(j,k,l)-un(j-1,k,l);
	
				pr(0)=pl(0);
				pr(1)=un(j+1,k,l)-un(j,k,l);
	
				pcx(0)=pl(0) - (1.0/12.0)*(un(j+1,k,l) + un(j,k+1,l) - 4*un(j,k,l) + un(j-1,k,l) + un(j,k-1,l));
				pcx(1)=.5*(pl(1)+pr(1));
				pcx(2)=pr(1)-pl(1);
				
				ISl_vector(l)=pl(1)*pl(1)/(u_norm(l)+eps);
				IScx_vector(l)=(1.0/(u_norm(l)+eps))*((13.0/3.0)*pcx(2)*pcx(2) + pcx(1)*pcx(1));
				ISr_vector(l)=pr(1)*pr(1)/(u_norm(l)+eps);
				
				pb(0)=pl(0);
				pb(1)=un(j,k,l)-un(j,k-1,l);
	
				pt(0)=pb(0);
				pt(1)=un(j,k+1,l)-un(j,k,l);
	
				pcy(0)=pcx(0);
				pcy(1)=.5*(pb(1)+pt(1));
				pcy(2)=pt(1)-pb(1);
				
				ISb_vector(l)=pb(1)*pb(1)/(u_norm(l)+eps);
				IScy_vector(l)=(1.0/(u_norm(l)+eps))*((13.0/3.0)*pcy(2)*pcy(2) + pcy(1)*pcy(1));
				ISt_vector(l)=pt(1)*pt(1)/(u_norm(l)+eps);
			}
			
			ISl(j,k)=(1.0/L)*sum_2d_SD3(ISl_vector);
			IScx(j,k)=(1.0/L)*sum_2d_SD3(IScx_vector);
			ISr(j,k)=(1.0/L)*sum_2d_SD3(ISr_vector);
			
			ISb(j,k)=(1.0/L)*sum_2d_SD3(ISb_vector);
			IScy(j,k)=(1.0/L)*sum_2d_SD3(IScy_vector);
			ISt(j,k)=(1.0/L)*sum_2d_SD3(ISt_vector);
		}
	}
}
