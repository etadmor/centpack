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

#include "centpack_2d_SD3.h"

using namespace std;

void CENTPACK::C_flux_2d_SD3(const doublearray3d& u_N, const doublearray3d& u_S, const doublearray3d& u_E, const doublearray3d& u_W, const doublearray3d& u_NE, const doublearray3d& u_SE, const doublearray3d& u_SW, const doublearray3d& u_NW, const double& lambda, const double& mu, const double& gamma, doublearray3d& C)
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
	doublearray1d u_nejm1(L);
	doublearray1d u_nekm1(L);
	doublearray1d u_ne(L);
	doublearray1d u_sejm1(L);
	doublearray1d u_se(L);
	doublearray1d u_sekp1(L);
	doublearray1d u_sw(L);
	doublearray1d u_swjp1(L);
	doublearray1d u_swkp1(L);
	doublearray1d u_nwkm1(L);
	doublearray1d u_nw(L);
	doublearray1d u_nwjp1(L);
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
				u_nejm1(l)=u_NE(j-1,k,l);
				u_nekm1(l)=u_NE(j,k-1,l);
				u_ne(l)=u_NE(j,k,l);
				u_sejm1(l)=u_SE(j-1,k,l);
				u_se(l)=u_SE(j,k,l);
				u_sekp1(l)=u_SE(j,k+1,l);
				u_sw(l)=u_SW(j,k,l);
				u_swjp1(l)=u_SW(j+1,k,l);
				u_swkp1(l)=u_SW(j,k+1,l);
				u_nwkm1(l)=u_NW(j,k-1,l);
				u_nw(l)=u_NW(j,k,l);
				u_nwjp1(l)=u_NW(j+1,k,l);
			}
		
			Hx_flux_2d_SD3(u_nw, u_w, u_sw, u_nejm1, u_ejm1, u_sejm1, Hx_halfm, gamma);
			Hx_flux_2d_SD3(u_nwjp1, u_wjp1, u_swjp1, u_ne, u_e, u_se, Hx_halfp, gamma);
			Hy_flux_2d_SD3(u_sw, u_s, u_se, u_nekm1, u_nkm1, u_nwkm1, Hy_halfm, gamma);
			Hy_flux_2d_SD3(u_swkp1, u_skp1, u_sekp1, u_ne, u_n, u_nw, Hy_halfp, gamma);
			
			for (l=0; l<L; l++)
				C(j,k,l)=-lambda*(Hx_halfp(l)-Hx_halfm(l)) - mu*(Hy_halfp(l)-Hy_halfm(l));
		}
	}
}
