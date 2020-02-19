////////////////////////////////////////////////////////////////////////////
//
// CentPack -- A generic numerical solver for hyperbolic conservation
//             laws and related time dependent problems
//
// Copyright (C) 2006 Jorge Balbas and Eitan Tadmor
//
// time_step_2d.cc -- a core function of CentPack
//
// Function called by centpack_2d_FD2() (main)
//
// Requires: spectral_radii
//
// This function calculates size of the time step to be used in the evolution of
// the solution variables at each time step
//
// Input (passed by reference):
// 
// (1) un  --  a doublearray3d type variable holding the cell averages of u
//              at time t=t^n over the discretized solution domain
// (2) dx  --  a double type variable holding the space scale dx
//
// (3) dy  --  a double type variable holding the space scale dy
//
// (4) cfl  --  a double type variable holding the CFL restriction for the 
//              central scheme implemented
//
// (5) gamma -- a double type variable holding the ratio of specific heats
//
// Output (returned by reference):
//
// (1) dt  --  a double type variable holding the value of the next time step to 
//             be used in next evolution step
//
// (2) lambda  --  a double type variable holding themesh ratio dt/dx
//
// (3) mu  --  a double type variable holding the mesh ratio dt/dy
//
// (4) t  --  a double type variable holding the simulation time completed so 
//            far
//
// (5) t_out  --  a double type variable holding the simulation time since last 
//                output
//
////////////////////////////////////////////////////////////////////////////////

#include "centpack_2d_FD2.h"
#include "centpack_2d_SD2.h"
#include "centpack_2d_SD3.h"

using namespace std;

void CENTPACK::time_step_2d(const doublearray3d& un, const double& dx, const double& dy, const double& cfl, double& dt, double& t, double& dt_out, double& t_out, double& lambda, double& mu, const double& gamma)
{
	long J=un.getIndex1Size() - 4;
	long K=un.getIndex2Size() - 4;
	long L=un.getIndex3Size();
	long j, k, l;
	double rx, ry, r_maxx, r_maxy;
	doublearray1d u_vector(L);
	
	r_maxx = 0.0;
	r_maxy = 0.0;
	
	for (j = 2; j < J+2; j++)
	{
		for (k = 2; k < K+2; k++)
		{
			for (l = 0; l < L; l++)
				u_vector(l) = un(j,k,l);
		
			spectral_radii(u_vector, gamma, rx, ry);
		
			if (rx > r_maxx)
				r_maxx = rx;
			if (ry > r_maxy)
				r_maxy = ry;
		}
	}
	
	dt = cfl/sqrt((r_maxx/dx)*(r_maxx/dx)+(r_maxy/dy)*(r_maxy/dy));
	dt = min(dt, dt_out - t_out);
	
	lambda = dt/dx;
	mu = dt/dy;
	
	t = t + dt;
	t_out = t_out + dt;
}
