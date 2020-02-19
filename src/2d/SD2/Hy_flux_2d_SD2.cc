////////////////////////////////////////////////////////////////////////////
//
// CentPack -- A generic numerical solver for hyperbolic conservation
//             laws and related time dependent problems
//
// Copyright (C) 2005 Jorge Balbas and Eitan Tadmor
//
// Hy_flux.cc -- a core function of CentPack
//
// Function called by C_flux
//
// Requires: flux_y and spectral_radii (see prototypes below)
//
// This function calculates the numerical flux of the semi-discre formulation 
// for the evolution of the cell averages of u at each side of the cell 
// interfaces y=y_k-1/2 and y=y_k+1/2
//
// Input (passed by reference):
//
// (1) two doublearray1d type variables
//
//     u_n  --  holds the point values of u_N for the cell centered at 
//              (x_j,z_k-1) or (x_j,z_k)
//
//     u_s  --  holds the point values of u_S for the cell centered at
//              (x_j,z_k) or (x_j,z_k+1)
//
// (2) gamma -- a double type variable holding the ratio of specific heats
//
// Output (returned by reference):
//
// (1) Hy -- numerical flux at cell interface y=y_k-1/2 or y=y_k+1/2
//
////////////////////////////////////////////////////////////////////////////////

#include "centpack_2d_SD2.h"

using namespace std;

void CENTPACK::Hy_flux_2d_SD2(const doublearray1d& u_s, const doublearray1d& u_n, doublearray1d& Hy, const double& gamma)
{
  long l;
  long L=Hy.getIndex1Size();
  double rw, re, rn, rs, b;

  doublearray1d g_n(L), g_s(L);
  
  spectral_radii(u_n, gamma, re, rn);
  spectral_radii(u_s, gamma, rw, rs);

  if (rn>rs)
	  b=rn;
  else
	  b=rs;
  
  flux_y(u_n, g_n, gamma);
  flux_y(u_s, g_s, gamma);
  
  for (l=0; l<L; l++)
	  Hy(l)=.5*((g_s(l)+g_n(l)) - b*(u_s(l)-u_n(l)));
}
