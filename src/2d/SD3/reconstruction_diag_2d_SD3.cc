////////////////////////////////////////////////////////////////////////////
//
// CentPack -- A generic numerical solver for hyperbolic conservation
//             laws and related time dependent problems
//
// Copyright (C) 2006 Jorge Balbas and Eitan Tadmor
//
// reconstruction_diag.cc -- a core function of CentPack
//
// Function called by evolution
//
// This function implements a 3rd order Kurganov-Levy CWENO reconstruction
// of the corner values of the solution u from its cell averages along
// the SW-NE and SE-NW directions
// 
// Input (passed by reference):
//
// (1) un  --  a doublearray3d type variable holding the cell averages of u
//             at time t=t^n over the discretized solution domain
//
// (2) Six doublearray2d arrays holding the value of the smoothness 
//     indicators used to calculate the non-linear weights of the 
//     non-oscillatory reconstruction along the diagonal directions
//     (i.e., SW-NE and SE-NW directions):
//
//     dISl  --  SW reconstruction
//
//     dIScx  --  central reconstruction along SW-NE direction
//
//     dISr  --  NE reconstruction
//
//     dISb  --  SE reconstruction
//
//     dIScy  --  central reconstruction along SE-NW direction
//
//     dISt  --  NW reconstruction
//
//
// Output (returned by reference):
//
// (1) Four doublearray3d type variables
//
//     u_NE  --  holds the point values of u at the cell interfaces
//               (x_j+1/2, y_k+1/2) over the entire solution domain
//
//     u_SE  --  holds the point values of u at the cell interfaces
//               (x_j+1/2, y_k-1/2) over the entire solution domain
//
//     u_SW  --  holds the point values of u at the cell interfaces
//               (x_j-1/2, y_k-1/2) over the entire solution domain
//
//     u_NW  --  holds the point values of u at the cell interfaces
//               (x_j-1/2, y_k+1/2) over the entire solution domain
//
////////////////////////////////////////////////////////////////////////////////

#include "centpack_2d_SD3.h"

using namespace std;

void CENTPACK::reconstruction_diag_2d_SD3(const doublearray3d& un, doublearray3d& u_NE, doublearray3d& u_SE, doublearray3d& u_SW, doublearray3d& u_NW, const doublearray2d& dISl, const doublearray2d& dIScx, const doublearray2d& dISr, const doublearray2d& dISb, const doublearray2d& dIScy, const doublearray2d& dISt)
{
	
	long J=un.getIndex1Size() - 4;
	long K=un.getIndex2Size() - 4;
	long L=un.getIndex3Size();
	long j, k, l, s;
	
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
				pl(1)=un(j,k,l)-un(j-1,k-1,l);
	
				pr(0)=pl(0);
				pr(1)=un(j+1,k+1,l)-un(j,k,l);
	
				pcx(0)=pl(0) - (1.0/12.0)*(un(j+1,k+1,l) + un(j-1,k+1,l) - 4*un(j,k,l) + un(j-1,k-1,l) + un(j+1,k-1,l));
				pcx(1)=.5*(pl(1)+pr(1));
				pcx(2)=pr(1)-pl(1);
				
				alpl=cl/((eps+dISl(j,k))*(eps+dISl(j,k)));
				alpcx=ccx/((eps+dIScx(j,k))*(eps+dIScx(j,k)));
				alpr=cr/((eps+dISr(j,k))*(eps+dISr(j,k)));
				
				alp_sum=alpl+alpcx+alpr;
				
				wl=alpl/alp_sum;
				wcx=alpcx/alp_sum;
				wr=alpr/alp_sum;
				
				pb(0)=pl(0);
				pb(1)=un(j,k,l)-un(j+1,k-1,l);
	
				pt(0)=pb(0);
				pt(1)=un(j-1,k+1,l)-un(j,k,l);
	
				pcy(0)=pcx(0);
				pcy(1)=.5*(pb(1)+pt(1));
				pcy(2)=pt(1)-pb(1);
				
				alpb=cb/((eps+dISb(j,k))*(eps+dISb(j,k)));
				alpcy=ccy/((eps+dIScy(j,k))*(eps+dIScy(j,k)));
				alpt=ct/((eps+dISt(j,k))*(eps+dISt(j,k)));
				
				alp_sum=alpb+alpcy+alpt;
				
				wb=alpb/alp_sum;
				wcy=alpcy/alp_sum;
				wt=alpt/alp_sum;
				
				u_NW(j,k,l)=wb*(pb(0)+.5*pb(1)) + wcy*(pcy(0)+.5*pcy(1)+.25*pcy(2)) +wt*(pt(0)+.5*pt(1));
				u_SE(j,k,l)=wb*(pb(0)-.5*pb(1)) + wcy*(pcy(0)-.5*pcy(1)+.25*pcy(2)) +wt*(pt(0)-.5*pt(1));
				u_NE(j,k,l)=wl*(pl(0)+.5*pl(1)) + wcx*(pcx(0)+.5*pcx(1)+.25*pcx(2)) +wr*(pr(0)+.5*pr(1));
				u_SW(j,k,l)=wl*(pl(0)-.5*pl(1)) + wcx*(pcx(0)-.5*pcx(1)+.25*pcx(2)) +wr*(pr(0)-.5*pr(1));
			}
		}
	}
}
