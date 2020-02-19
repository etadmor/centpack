////////////////////////////////////////////////////////////////////////////
//
// CentPack -- A generic numerical solver for hyperbolic conservation
//             laws and related time dependent problems
//
// Copyright (C) 2006 Jorge Balbas and Eitan Tadmor
//
// C_flux.cc -- a core function of CentPack
//
// Function called by evolution
//
// Requires: Hx_flux and Hy_flux (see prototypes below)
//
// This function calculates the SSP-RK fluxes employed for the evolution of the
// cell averages of u
//
// Input (passed by reference):
//
// (1) Eight doublearray3d type variables:
//
//     u_N  --  holds the point values of u at the cell interfaces
//                (x_j, y_k+1/2) over the entire solution domain
//
//     u_S  --  holds the point values of u at the cell interfaces
//                (x_j, y_k-1/2) over the entire solution domain
//
//     u_E  --  holds the point values of u at the cell interfaces
//                (x_j+1/2, y_k) over the entire solution domain
//
//     u_W  --  holds the point values of u at the cell interfaces
//                (x_j-1/2, y_k) over the entire solution domain
//
// (2) lambda  --  a double type variable holding the mesh ratio dt/dx
//
// (3) mu  --  a double type variable holding the mesh ratio dt/dy
//
// (4) gamma -- a double type variable holding the ratio of specific heats
//
// Output (returned by reference):
//
// (1) C -- a doublearray3d type variable holding the value of the SSP-RK fluxes
//
////////////////////////////////////////////////////////////////////////////////

#include "centpack_2d_SD2.h"

using namespace std;

void CENTPACK::C_flux_2d_SD2(const doublearray3d& u_N, const doublearray3d& u_S, const doublearray3d& u_E, const doublearray3d& u_W, const double& lambda, const double& mu, const double& gamma, doublearray3d& C)
{
	
	long j, k, l, J, K, L;
	
	J=C.getIndex1Size() - 4;
	K=C.getIndex2Size() - 4;
	L=C.getIndex3Size();
	
	doublearray1d u_nkm1(L);
	doublearray1d u_n(L);
	doublearray1d u_s(L);
	doublearray1d u_skp1(L);
	doublearray1d u_ejm1(L);
	doublearray1d u_e(L);
	doublearray1d u_w(L);
	doublearray1d u_wjp1(L);
	
	doublearray1d Hx_halfp(L), Hx_halfm(L), Hy_halfp(L), Hy_halfm(L);

	for (k=2; k<K+2; k++)
	{
		for (j=2; j<J+2; j++)
		{
			for (l=0; l<L; l++)
			{
				u_nkm1(l)=u_N(j,k-1,l);
				u_n(l)=u_N(j,k,l);
				u_s(l)=u_S(j,k,l);
				u_skp1(l)=u_S(j,k+1,l);
				u_ejm1(l)=u_E(j-1,k,l);
				u_e(l)=u_E(j,k,l);
				u_w(l)=u_W(j,k,l);
				u_wjp1(l)=u_W(j+1,k,l);
			}
		
			Hx_flux_2d_SD2(u_w, u_ejm1, Hx_halfm, gamma);
			Hx_flux_2d_SD2(u_wjp1, u_e, Hx_halfp, gamma);
			Hy_flux_2d_SD2(u_s, u_nkm1, Hy_halfm, gamma);
			Hy_flux_2d_SD2(u_skp1, u_n, Hy_halfp, gamma);
			
			for (l=0; l<L; l++)
				C(j,k,l)=-lambda*(Hx_halfp(l)-Hx_halfm(l)) - mu*(Hy_halfp(l)-Hy_halfm(l));
		}
	}
}
