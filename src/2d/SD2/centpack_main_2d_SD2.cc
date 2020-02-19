////////////////////////////////////////////////////////////////////////////
//
// CentPack -- A generic numerical solver for hyperbolic conservation
//             laws and related time dependent problems
//
// Copyright (C) 2006 Jorge Balbas and Eitan Tadmor
//
// solver.cc -- main function of CentPack
//
// Requires: read_in, initial_conditions.cc, evolution.cc, writeout.cc
// and those called by these (see README file for more information)f
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
//               y-interval of the solution
//
// (4) y_final -- a double type variable holding the right end point of the 
//                y-interval of the solution
//
// (5) J  --  a long type variable holding the number of cells along x-dimension
//
// (6) K  --  a long type variable holding the number of cells along y-dimension
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
// (11) dt_out  --  a double type variable holding the desired output interval
//
// (12) alpha  --  minmod limiting parameter
//
// Output:
//
// This function does not produce any direct output except for monitoring 
// information (see README file).  The solutions variables are output via a 
// writeout function called at the time interval indicated by the user in the 
// input file
//
////////////////////////////////////////////////////////////////////////////////

#include "centpack_2d_SD2.h"

using namespace std;

int CENTPACK::centpack_main_2d_SD2()
{
	
	double x_init, x_final, y_init, y_final, t, t_init, t_final, t_out;
	double dx, dy, dt, cfl, dt_out, alpha;
	double rx, ry, r_maxx, r_maxy, lambda, mu;
	double gamma;
	double sum_t=0.0;
	double dt_cpu=0.0;
	double t_start=clock();
	double pi=4.0*atan2(1.0,1.0);
	long J, K, L;
	long j, k, l;
	long n=0;
	
	read_in_2d(x_init, x_final, y_init, y_final, J, K, L, gamma, t_final, dt_out, cfl, alpha);
	
	doublearray1d x(J+4), y(K+4);
	
	cout.setf(ios::scientific, ios::floatfield);
	
	t=0.0;
	t_init=0.0;
	dx=(x_final - x_init)/J;
	dy=(y_final-y_init)/K;
	
	doublearray1d u_vector(L);
	doublearray3d un(J+4,K+4,L);
	
	initial_conditions(un, x_init, x_final, y_init, y_final, dx, dy, gamma, x, y);
	
	writeout(un,t,gamma,n);
	t_out=0.0;
	n++;
	
	do
	{
		time_step_2d(un, dx, dy, cfl, dt, t, dt_out, t_out, lambda, mu, gamma);
		evolution_2d_SD2(un, lambda, mu, dx, dy, gamma, alpha);
		
		dt_cpu=(clock()-t_start)/CLOCKS_PER_SEC;
		sum_t=sum_t+dt_cpu;
		
		end_of_step_2d_SD2(un, dt, t, dt_out, t_out, n, gamma, sum_t, dt_cpu);
		
		t_start=clock();
		
	}while(t<t_final);
	
	writeout(un,t,gamma,n);
	
	run_info_2d(dt, sum_t, J, K, cfl);
	
	return 0;

}
//////////////////////////////////////////////////////////////////////////////
