////////////////////////////////////////////////////////////////////////////
//
// CentPack -- A generic numerical solver for hyperbolic conservation
//             laws and related time dependent problems
//
// Copyright (C) 2006 Jorge Balbas and Eitan Tadmor
//
// reconstruction.cc -- a core function of CentPack
//
// Function called by evolution
//
// This function implements a 3rd order Kurganov-Levy CWENO reconstruction
// of the interface values of the solution u from its cell averages along
// the x- and y-directions
// 
// Input (passed by reference):
//
// (1) un  --  a doublearray3d type variable holding the cell averages of u
//             at time t=t^n over the discretized solution domain
//
// (2) Six doublearray2d arrays holding the value of the smoothness 
//     indicators used to calculate the non-linear weights of the 
//     non-oscillatory reconstruction along the x- and y-directions:
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

#include "centpack_2d_SD3.h"

using namespace std;

void CENTPACK::reconstruction_2d_SD3(doublearray3d& un, doublearray3d& u_N, doublearray3d& u_S, doublearray3d& u_E, doublearray3d& u_W, const doublearray2d& ISl, const doublearray2d& IScx, const doublearray2d& ISr, const doublearray2d& ISb, const doublearray2d& IScy, const doublearray2d& ISt)
{
	
	long J=un.getIndex1Size() - 4;
	long K=un.getIndex2Size() - 4;
	long L=un.getIndex3Size();
	long j, k, l;
	
	doublearray1d pl(2), pcx(3), pr(2);
	doublearray1d pt(2), pcy(3), pb(2);
	double wl, wcx, wr, wt, wcy, wb; 
	double alpl, alpcx, alpr, alpb, alpcy, alpt, alp_sum;
	double cl=.25, ccx=.5, cr=.25;
	double cb=.25, ccy=.5, ct=.25;
	double eps=0.000001;
	
	for (l=0; l<L; l++)
	{		
		for (k=2; k<K+2; k++)
		{
			for (j=2; j<J+2; j++)
			{
				pl(0)=un(j,k,l);
				pl(1)=un(j,k,l)-un(j-1,k,l);
	
				pr(0)=pl(0);
				pr(1)=un(j+1,k,l)-un(j,k,l);
	
				pcx(0)=pl(0) - (1.0/12.0)*(un(j+1,k,l) + un(j,k+1,l) - 4*un(j,k,l) + un(j-1,k,l) + un(j,k-1,l));
				pcx(1)=.5*(pl(1)+pr(1));
				pcx(2)=pr(1)-pl(1);
				
				alpl=cl/((eps+ISl(j,k))*(eps+ISl(j,k)));
				alpcx=ccx/((eps+IScx(j,k))*(eps+IScx(j,k)));
				alpr=cr/((eps+ISr(j,k))*(eps+ISr(j,k)));
				
				alp_sum=alpl+alpcx+alpr;
				
				wl=alpl/alp_sum;
				wcx=alpcx/alp_sum;
				wr=alpr/alp_sum;
				
				pb(0)=pl(0);
				pb(1)=un(j,k,l)-un(j,k-1,l);
	
				pt(0)=pb(0);
				pt(1)=un(j,k+1,l)-un(j,k,l);
	
				pcy(0)=pcx(0);
				pcy(1)=.5*(pb(1)+pt(1));
				pcy(2)=pt(1)-pb(1);
				
				alpb=cb/((eps+ISb(j,k))*(eps+ISb(j,k)));
				alpcy=ccy/((eps+IScy(j,k))*(eps+IScy(j,k)));
				alpt=ct/((eps+ISt(j,k))*(eps+ISt(j,k)));
				
				alp_sum=alpb+alpcy+alpt;
				
				wb=alpb/alp_sum;
				wcy=alpcy/alp_sum;
				wt=alpt/alp_sum;
				
				u_N(j,k,l)=wb*(pb(0)+.5*pb(1)) + wcy*(pcy(0)+.5*pcy(1)+.25*pcy(2)) +wt*(pt(0)+.5*pt(1));
				u_S(j,k,l)=wb*(pb(0)-.5*pb(1)) + wcy*(pcy(0)-.5*pcy(1)+.25*pcy(2)) +wt*(pt(0)-.5*pt(1));
				u_E(j,k,l)=wl*(pl(0)+.5*pl(1)) + wcx*(pcx(0)+.5*pcx(1)+.25*pcx(2)) +wr*(pr(0)+.5*pr(1));
				u_W(j,k,l)=wl*(pl(0)-.5*pl(1)) + wcx*(pcx(0)-.5*pcx(1)+.25*pcx(2)) +wr*(pr(0)-.5*pr(1));
			}
		}
	}
}
