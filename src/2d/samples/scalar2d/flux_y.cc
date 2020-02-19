//////////////////////////////////////////////////////////////////////////
//
// CentPack -- A generic numerical solver for hyperbolic conservation
//             laws and related time dependent problems
//
// Copyright (C) 2006 Jorge Balbas and Eitan Tadmor
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
	g(0) = (1.0/3.0)*pow(u(0),3.0);
}
