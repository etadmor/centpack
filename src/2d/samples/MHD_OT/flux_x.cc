////////////////////////////////////////////////////////////////////////////////
//
// CentPack -- A generic numerical solver for hyperbolic conservation
//             laws and related time dependent problems
//
// Copyright (C) 2006 Jorge Balbas and Eitan Tadmor
// 
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the warranty of MERCHANTABILITY or FITNESS FOR A
// PARTICULAR PURPOSE.
//
// flux_x.cc -- Ideal MHD x-flux function, an auxiliary function of
//
// Function called by Hx_flux, (model specific)
//
// Input (passed by reference):
//
// (1) u -- a doublearray1d type variable holding the conserved
//          quantities density, momentum (3 components), magnetic field
//          (3 components), and total energy at a point of the didcretized 
//          solution domain.
//
// (2) gamma -- a double type variable holding the ratio of specific heats
//
// Output (returned by reference):
//
// (1) f -- a doublearray1d type variable holding f(u), the flux in
//          the x direction.
//
////////////////////////////////////////////////////////////////////////////////

#include"centpack_2d_FD2.h"
#include"centpack_2d_SD2.h"

using namespace std;

void CENTPACK::flux_x(const doublearray1d& u, doublearray1d& f, const double& gamma)
{
  double p, p_star;

  p=(gamma-1.0)*(u(7)-.5*((pow(u(1),2)+pow(u(2),2)+pow(u(3),2))/u(0)) - .5*(pow(u(4),2)+pow(u(5),2)+pow(u(6),2)));
  p_star= p + .5*(pow(u(4),2)+pow(u(5),2)+pow(u(6),2));
  
  f(0)=u(1);
  f(1)=(pow(u(1),2))/u(0) + p_star - pow(u(4),2);
  f(2)=(u(1)*u(2)/u(0)) - u(4)*u(5);
  f(3)=(u(1)*u(3)/u(0)) - u(4)*u(6);
  f(4)=0.0;
  f(5)=(u(5)*u(1)/u(0)) - u(4)*u(2)/u(0);
  f(6)=(u(6)*u(1)/u(0)) - u(4)*u(3)/u(0);
  f(7)=(u(7)+p_star)*(u(1)/u(0)) - u(4)*(u(4)*(u(1)/u(0)) + (u(2)/u(0))*u(5)+(u(3)/u(0))*u(6));
}
