////////////////////////////////////////////////////////////////////////////////
//
// CentPack -- A generic numerical solver for hyperbolic conservation
//             laws and related time dependent problems
//
// Copyright (C) 2006 Jorge Balbas and Eitan Tadmor
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
#include"centpack_2d_SD3.h"

using namespace std;

void CENTPACK::flux_x(const doublearray1d& u, doublearray1d& f, const double& gamma)
{
	f(0) = sin(u(0));
}
