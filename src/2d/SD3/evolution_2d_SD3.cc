////////////////////////////////////////////////////////////////////////////
//
// CentPack -- A generic numerical solver for hyperbolic conservation
//             laws and related time dependent problems
//
// Copyright (C) 2006 Jorge Balbas and Eitan Tadmor
//
// evolution.cc -- a core function of CentPack
//
// Function called by solver (main)
//
// Requires: boundary_conditions, norm, reconstruction, reconstruction_diag,
// indicators, indicators_diag and C_flux (see prototypes below)
//
// This function evolves the cell averages of u from t=t^n to t=t^n + dt
// 
// Input (passed by reference):
//
// (1)  un  --  a doublearray3d type variable holding the cell averages of u
//              at time t=t^n over the discretized solution domain
//
// (2) dx  --  a double type variable holding the space scale dx
//
// (3) dy  --  a double type variable holding the space scale dy
//
// (4) lambda  --  a double type variable holding the mesh ratio dt/dx
//
// (5) mu  --  a double type variable holding the mesh ratio dt/dy
//
// (6) gamma -- a double type variable holding the ratio of specific heats
//
// Output (returned by reference):
//
// (1)  un  --  a doublearray3d type variable holding the cell averages of u
//              at time t=t^n + dt over the discretized solution domain
//
////////////////////////////////////////////////////////////////////////////////

#include "centpack_2d_SD3.h"

using namespace std;

void CENTPACK::evolution_2d_SD3(doublearray3d& un, const double& lambda, const double& mu, const double& dx, const double& dy, const double& gamma)
{
	long j, k, l, J, K, L;
	
	J=un.getIndex1Size() - 4;
	K=un.getIndex2Size() - 4;
	L=un.getIndex3Size();
	
	doublearray1d u_norm(L);
	
	doublearray2d ISl(J+4,K+4), IScx(J+4,K+4), ISr(J+4,K+4), ISb(J+4,K+4), IScy(J+4,K+4), ISt(J+4,K+4);
	
	doublearray2d dISl(J+4,K+4), dIScx(J+4,K+4), dISr(J+4,K+4), dISb(J+4,K+4), dIScy(J+4,K+4), dISt(J+4,K+4);
	
	doublearray3d u(J+4,K+4,L);
	doublearray3d u_N(J+4,K+4,L);
	doublearray3d u_S(J+4,K+4,L);
	doublearray3d u_E(J+4,K+4,L);
	doublearray3d u_W(J+4,K+4,L);
	doublearray3d u_NE(J+4,K+4,L);
	doublearray3d u_SE(J+4,K+4,L);
	doublearray3d u_SW(J+4,K+4,L);
	doublearray3d u_NW(J+4,K+4,L);
	doublearray3d C0(J+4,K+4,L), C1(J+4,K+4,L), C2(J+4,K+4,L);
	
	norm_2d_SD3(un, u_norm, dx, dy);
	boundary_conditions(un);
	indicators_2d_SD3(un, ISl, IScx, ISr, ISb, IScy, ISt, u_norm);
	reconstruction_2d_SD3(un, u_N, u_S, u_E, u_W, ISl, IScx, ISr, ISb, IScy, ISt);
	indicators_diag_2d_SD3(un, dISl, dIScx, dISr, dISb, dIScy, dISt, u_norm);
	reconstruction_diag_2d_SD3(un, u_NE, u_SE, u_SW, u_NW, dISl, dIScx, dISr, dISb, dIScy, dISt);
	boundary_conditions(u_N);
	boundary_conditions(u_S);
	boundary_conditions(u_E);
	boundary_conditions(u_W);
	boundary_conditions(u_NE);
	boundary_conditions(u_SE);
	boundary_conditions(u_SW);
	boundary_conditions(u_NW);
	C_flux_2d_SD3(u_N, u_S, u_E, u_W, u_NE, u_SE, u_SW, u_NW, lambda, mu, gamma, C0);
	
	for (l=0; l<L; l++)
	{
		for (k=2; k<K+2; k++)
		{
			for (j=2; j<J+2; j++)
				un(j,k,l)=un(j,k,l)+C0(j,k,l);
		}
	}
	
	reconstruction_2d_SD3(un, u_N, u_S, u_E, u_W, ISl, IScx, ISr, ISb, IScy, ISt);
	reconstruction_diag_2d_SD3(un, u_NE, u_SE, u_SW, u_NW, dISl, dIScx, dISr, dISb, dIScy, dISt);
	boundary_conditions(u_N);
	boundary_conditions(u_S);
	boundary_conditions(u_E);
	boundary_conditions(u_W);
	boundary_conditions(u_NE);
	boundary_conditions(u_SE);
	boundary_conditions(u_SW);
	boundary_conditions(u_NW);
	C_flux_2d_SD3(u_N, u_S, u_E, u_W, u_NE, u_SE, u_SW, u_NW, lambda, mu, gamma, C1);
	
	for (l=0; l<L; l++)
	{
		for (k=2; k<K+2; k++)
		{
			for (j=2; j<J+2; j++)
				un(j,k,l)=un(j,k,l)+.25*(-3*C0(j,k,l)+C1(j,k,l));
		}
	}
	
	reconstruction_2d_SD3(un, u_N, u_S, u_E, u_W, ISl, IScx, ISr, ISb, IScy, ISt);
	reconstruction_diag_2d_SD3(un, u_NE, u_SE, u_SW, u_NW, dISl, dIScx, dISr, dISb, dIScy, dISt);
	boundary_conditions(u_N);
	boundary_conditions(u_S);
	boundary_conditions(u_E);
	boundary_conditions(u_W);
	boundary_conditions(u_NE);
	boundary_conditions(u_SE);
	boundary_conditions(u_SW);
	boundary_conditions(u_NW);
	C_flux_2d_SD3(u_N, u_S, u_E, u_W, u_NE, u_SE, u_SW, u_NW, lambda, mu, gamma, C2);
	
	for (l=0; l<L; l++)
	{
		for (k=2; k<K+2; k++)
		{
			for (j=2; j<J+2; j++)
				un(j,k,l)=un(j,k,l) + (1.0/12.0)*(-C0(j,k,l) - C1(j,k,l) + 8*C2(j,k,l));
		}
	}
	
	boundary_conditions(un);
}
