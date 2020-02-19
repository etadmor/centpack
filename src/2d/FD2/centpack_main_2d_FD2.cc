////////////////////////////////////////////////////////////////////////////
//
// CentPack -- A generic numerical solver for hyperbolic conservation
//             laws and related time dependent problems
//
// Copyright (C) 2006 Jorge Balbas and Eitan Tadmor
// 
// centpack_main_2d_FD2.cc -- main function in
//
// Requires: initial_conditions.cc, spectral_radii.cc, evolution.cc, writeout.cc
// and those called by these (see README file for more information.
//
// This program approximates the solution of hyperbolic conservation laws in 
// two space dimensions (denoted by x and z),
//
//   u_t + f(u)_x + g(u)_y = 0                 (1)
//
// Using Jiang and Tadmor's 2nd order central scheme
//
// Input (read from file input):
// 
// (1) x_init -- a double type variable holding the left end point of the 
//               x-interval of the solution
//
// (2) x_final -- a double type variable holding the right end point of the 
//                x-interval of the solution
//
// (3) z_init -- a double type variable holding the left end point of the 
//               z-interval of the solution
//
// (4) z_final -- a double type variable holding the right end point of the 
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
// (12) alpha  --  a double type variable holding the minmod3 parameter
//
// Output:
//
// This function does not produce any direct output except for monitoring 
// information (see README file).  The solutions variables are output via a 
// writeout function called at the time interval indicated by the user in the 
// input file
//
////////////////////////////////////////////////////////////////////////////////

#include "centpack_2d_FD2.h"

using namespace std;

int CENTPACK::centpack_main_2d_FD2()
{
	disclaimer();
	
	double x_init, x_final, y_init, y_final, t_final;
	double dx, dy, dt, cfl, dt_out;
	double lambda, mu;
	double gamma, alpha;
	double sum_t=0.0;
	double dt_cpu=0.0;
	double t_start=clock();
	double t=0.0;
	double t_out=0.0;
	long J, K, L;
	long n=0;
	bool odd=true;
	
	read_in_2d(x_init, x_final, y_init, y_final, J, K, L, gamma, t_final, dt_out, cfl, alpha);
	
	cout.setf(ios::scientific, ios::floatfield);
	
	doublearray1d x(J+4), y(K+4);
	doublearray3d un(J+4,K+4,L), unhalf(J+4,K+4,L);
	
	initial_conditions(un, x_init, x_final, y_init, y_final, dx, dy, gamma, x, y);
	writeout(un, t, gamma, n);
	n++;
	
	do
	{
		time_step_2d(un, dx, dy, cfl, dt, t, dt_out, t_out, lambda, mu, gamma);
		reconstruction_2d_FD2(un, unhalf, alpha, odd);
		evolution_2d_FD2(un, unhalf, lambda, mu, dx, dy, gamma, alpha, odd);
		
		dt_cpu=(clock()-t_start)/CLOCKS_PER_SEC;
		sum_t=sum_t+dt_cpu;
		
		end_of_step_2d_FD2(un, dt, t, dt_out, t_out, n, gamma, sum_t, dt_cpu, odd);
		
		t_start=clock();
		
	}while(t<t_final);
	
	writeout(un, t, gamma, n);
	run_info_2d(dt, sum_t, J, K, cfl);
	
	return 0;

}
//////////////////////////////////////////////////////////////////////////////
