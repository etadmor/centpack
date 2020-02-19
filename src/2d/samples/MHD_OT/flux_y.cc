//////////////////////////////////////////////////////////////////////////
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
// flux_y.cc -- Ideal MHD y-flux function, an auxiliary function of
//
// Function called by Hy_flux
//
// Input (passed by reference):
//
// (1) u -- a doublearray1d type variable holding the conserved
//          quantities density, momentum (3 components), magnetic field
//          (3 components), and total energy at a point of the discretized 
//          solution domain.
//
// (2) gamma -- a double type variable holding the ratio of specific heats
//
// Output (returned by reference):
//
// (1) g -- a doublearray1d type variable holding g(u), the flux in
//          the y direction.
//
//////////////////////////////////////////////////////////////////////////

#include"centpack_2d_FD2.h"
#include"centpack_2d_SD2.h"
#include"centpack_2d_SD3.h"

using namespace std;

void CENTPACK::flux_y(const doublearray1d& u, doublearray1d& g, const double& gamma)
{
  double p, p_star;

  p=(gamma-1.0)*(u(7)-.5*((pow(u(1),2)+pow(u(2),2)+pow(u(3),2))/u(0)) - .5*(pow(u(4),2)+pow(u(5),2)+pow(u(6),2)));
  p_star= p + .5*(pow(u(4),2)+pow(u(5),2)+pow(u(6),2));
  
  g(0)=u(3);
  g(1)=(u(1)*u(3)/u(0)) - u(4)*u(6);
  g(2)=(u(2)*u(3)/u(0)) - u(5)*u(6);
  g(3)=(pow(u(3),2))/u(0) + p_star - pow(u(6),2);
  g(4)=(u(4)*u(3)/u(0)) - u(6)*u(1)/u(0);
  g(5)=(u(5)*u(3)/u(0)) - u(6)*u(2)/u(0);
  g(6)=0.0;
  g(7)=(u(7)+p_star)*(u(3)/u(0)) - u(6)*(u(4)*(u(1)/u(0)) + (u(2)/u(0))*u(5)+(u(3)/u(0))*u(6));
}
