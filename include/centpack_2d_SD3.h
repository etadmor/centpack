////////////////////////////////////////////////////////////////////////////
//
// CentPack -- A generic numerical solver for hyperbolic conservation
//             laws and related time dependent problems
//
// Copyright (C) 2006 Jorge Balbas and Eitan Tadmor
//
// centpack_main_2d_SD3.cc -- main function of CentPack
//
// Requires: read_in, initial_conditions, evolution.cc, writeout.cc
// and those called by these (see README file for more information).
//
// This program approximates the solution of hyperbolic conservation laws in 
// two space dimensions (denoted by x and y),
//
//   u_t + f(u)_x + g(u)_y = 0                 (1)
//
// Input (read from file input):
// 
// (1) x_init -- a double type variable holding the left end point of the 
//               x-interval of the solution
//
// (2) x_final -- a double type variable holding the right end point of the 
//                x-interval of the solution
//
// (3) y_init -- a double type variable holding the left end point of the 
//               z-interval of the solution
//
// (4) y_final -- a double type variable holding the right end point of the 
//                z-interval of the solution
//
// (5) J  --  a long type variable holding the number of cells along x-dimension
//
// (6) K  --  a long type variable holding the number of cells along z-dimension
//
// (7) L  --  a long type variable holding the number of components in the 
//            system
//
// (8) gamma - a double type variable holding the ratio of specific heats
//
// (9) cfl  --  a double type variable holding the CFL restriction to determine 
//              time step
// 
// (10) t_final  --  a double type variable holding the time of simulation
//
// (11) dt_frame  --  a double type variable holding the desired output interval
//
// Output:
//
// This function does not produce any direct output except for monitoring 
// information (see README file).  The solutions variables are output via a 
// writeout function called at the time interval indicated by the user in the 
// input file
//
////////////////////////////////////////////////////////////////////////////////

#ifndef CENTPACK_2d_SD3_H
#define CENTPACK_2d_SD3_H

#include "arrays.h"
#include<algorithm>
#include<cmath>
#include<cstdio>
#include<cstdlib>
#include<ctime>
#include<fstream>
#include<string>
#include<iostream>

namespace CENTPACK
{

////////////////////////////////////////////////////////////////////////////////
// MAIN CENTPACK FUNCTION 
////////////////////////////////////////////////////////////////////////////////
	
	int centpack_main_2d_SD3();
	
	void disclaimer();
	
	void read_in_2d(double& x_init, double& x_final, double& y_init, double& y_final, long& J, long& K, long& L, double& gamma, double& t_final, double& dt_out, double& cfl, double& alpha);
	
	void time_step_2d(const doublearray3d& un, const double& dx, const double& dy, const double& cfl, double& dt, double& t, double& dt_out, double& t_out, double& lambda, double& mu, const double& gamma);
	
	////////////////////////////////////////////////////////////////////////////
	// EVOLUTION
	////////////////////////////////////////////////////////////////////////////

	void evolution_2d_SD3(doublearray3d& un, const double& lambda, const double& mu, const double& dx, const double& dy, const double& gamma);
	
	void norm_2d_SD3(const doublearray3d& un, doublearray1d& u_norm, const double& dx, const double& dy);
	
	void indicators_2d_SD3(doublearray3d& un, doublearray2d& ISl, doublearray2d& IScx, doublearray2d& ISr, doublearray2d& ISb, doublearray2d& IScy, doublearray2d& ISt, const doublearray1d& u_norm);
	
	void indicators_diag_2d_SD3(const doublearray3d& un, doublearray2d& dISl, doublearray2d& dIScx, doublearray2d& dISr, doublearray2d& dISb, doublearray2d& dIScy, doublearray2d& dISt, const doublearray1d& u_norm);
	
	double sum_2d_SD3(const doublearray1d& u);

	void reconstruction_2d_SD3(doublearray3d& un, doublearray3d& u_N, doublearray3d& u_S, doublearray3d& u_E, doublearray3d& u_W, const doublearray2d& ISl, const doublearray2d& IScx, const doublearray2d& ISr, const doublearray2d& ISb, const doublearray2d& IScy, const doublearray2d& ISt);
	
	void reconstruction_diag_2d_SD3(const doublearray3d& un, doublearray3d& u_NE, doublearray3d& u_SE, doublearray3d& u_SW, doublearray3d& u_NW, const doublearray2d& dISl, const doublearray2d& dIScx, const doublearray2d& dISr, const doublearray2d& dISb, const doublearray2d& dIScy, const doublearray2d& dISt);
	
	void C_flux_2d_SD3(const doublearray3d& u_N, const doublearray3d& u_S, const doublearray3d& u_E, const doublearray3d& u_W, const doublearray3d& u_NE, const doublearray3d& u_SE, const doublearray3d& u_SW, const doublearray3d& u_NW, const double& lambda, const double& mu, const double& gamma, doublearray3d& C);
	
	void Hx_flux_2d_SD3(const doublearray1d& u_nw, const doublearray1d& u_w, const doublearray1d& u_sw, const doublearray1d& u_ne, const doublearray1d& u_e, const doublearray1d& u_se, doublearray1d& Hx, const double& gamma);
	
	void Hy_flux_2d_SD3(const doublearray1d& u_sw, const doublearray1d& u_s, const doublearray1d& u_se, const doublearray1d& u_ne, const doublearray1d& u_n, const doublearray1d& u_nw, doublearray1d& Hy, const double& gamma);
	
	//end of evolution routines

	void end_of_step_2d_SD3(const doublearray3d& un, const double& dt, const double& t, double& dt_out, double& t_out, long& n, const double& gamma, const double& sum_t, const double& dt_cpu);
	
	void run_info_2d(double& dt, double& sum_t, long& J, long& K, double& cfl);
	
////////////////////////////////////////////////////////////////////////////////
// AUXILIARY FUNCTIONS -- User defined, see examples
////////////////////////////////////////////////////////////////////////////////

	void boundary_conditions(doublearray3d& u);

	void flux_x(const doublearray1d& u, doublearray1d& f, const double& gamma);

	void flux_y(const doublearray1d& u, doublearray1d& g, const double& gamma);
	
	void initial_conditions(doublearray3d& un, const double& x_init, const double& x_final, const double& y_init, const double& y_final, double& dx, double& dy, const double& gamma, doublearray1d& x, doublearray1d& y);

	void spectral_radii(const doublearray1d& u, const double& gamma, double& rx, double& ry);

	void writeout(const doublearray3d& un, const double& t, const double& gamma, const long& n);
}

#endif
