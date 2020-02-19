////////////////////////////////////////////////////////////////////////////
//
// CentPack -- A generic numerical solver for hyperbolic conservation
//             laws and related time dependent problems
//
// Copyright (C) 2006 Jorge Balbas and Eitan Tadmor
//
// Hx_flux.cc -- a core function of CentPack
//
// Function called by C_flux
//
// Requires: flux_x and spectral_radii (see prototypes below)
//
// This function calculates the numerical flux of the semi-discre formulation 
// for the evolution of the cell averages of u at each side of the cell 
// interfaces x=x_j-1/2 and x=x_j+1/2
//
// Input (passed by reference):
//
// (1) six doublearray1d type variables
//
//     u_e  --  holds the point values of u_W for the cell centered at 
//              (x_j-1,y_k) or (x_j,y_k)
//
//     u_w  --  holds the point values of u_S for the cell centered at
//              (x_j,y_k) or (x_j+1,y_k)
//
//     u_ne  --  holds the point values of u_NE for the cell centered at
//               (x_j,y_k-1) or (x_j,y_k)
//
//     u_se  --  holds the point values of u_SE for the cell centered at
//               (x_j,y_k) or (x_j, y_k+1)
//
//     u_sw  --  holds the point values of u_SW for the cell centered at
//               (x_j,y_k) or (x_j,y_k+1)
//
//     u_nw  --  holds the point values of u_NW for the cell centered at
//               (x_j,y_k-1) or (x_j,y_k)
//
// (2) gamma -- a double type variable holding the ratio of specific heats
//
// Output (returned by reference):
//
// (1) Hx -- numerical flux at cell interface x=x_j-1/2 or x=x_j+1/2
//
////////////////////////////////////////////////////////////////////////////////

#include "centpack_2d_SD3.h"

using namespace std;

void CENTPACK::Hx_flux_2d_SD3(const doublearray1d& u_nw, const doublearray1d& u_w, const doublearray1d& u_sw, const doublearray1d& u_ne, const doublearray1d& u_e, const doublearray1d& u_se, doublearray1d& Hx, const double& gamma)
{
  long l;
  long L=Hx.getIndex1Size();
  double rw, re, rn, rs, a;
  
  doublearray1d f_e(L), f_w(L), f_ne(L), f_nw(L), f_se(L), f_sw(L);
  
  spectral_radii(u_e, gamma, re, rn);
  spectral_radii(u_w, gamma, rw, rs);

  if (re>rw)
	  a=re;
  else
	  a=rw;
  
  flux_x(u_e, f_e, gamma);
  flux_x(u_w, f_w, gamma);
  flux_x(u_ne, f_ne, gamma);
  flux_x(u_nw, f_nw, gamma);
  flux_x(u_se, f_se, gamma);
  flux_x(u_sw, f_sw, gamma);
  
  for (l=0; l<L; l++)
	  Hx(l)=(1.0/12.0)*((f_nw(l)+f_ne(l)+4.0*(f_w(l)+f_e(l))+f_sw(l)+f_se(l)) - a*(u_nw(l)-u_ne(l)+4.0*(u_w(l)-u_e(l))+u_sw(l)-u_se(l)));
}
