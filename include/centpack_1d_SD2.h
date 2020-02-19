#ifndef CENTPACK_1d_SD2_H
#define CENTPACK_1d_SD2_H

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
	
	int centpack_main_1d_SD2();

	void disclaimer();
	
	void read_in_1d(double& x_init, double& x_final, long& J, long& L, double& gamma, double& t_final, double& dt_out, double& cfl, double& alpha, double& B1);
	
	void time_step_1d(const doublearray2d& un, const double& dx, const double& cfl, double& dt, double& t, double& dt_out, double& t_out, double& lambda, const double& gamma, const double& B1);
	
	////////////////////////////////////////////////////////////////////////////
	// EVOLUTION
	////////////////////////////////////////////////////////////////////////////
	
	void evolution_1d_SD2(doublearray2d& un, const double& lambda, const double& dx, const double& gamma, const double& alpha, const double& B1);
	
	// for SD formulation, reconstruction is embeded in evolution
	
	void reconstruction_1d_SD2(doublearray2d& un, doublearray2d& u_E, doublearray2d& u_W, const double& alpha);
	
	void C_flux_1d_SD2(const doublearray2d& u_E, const doublearray2d& u_W, const double& lambda, const double& gamma, const double& B1, doublearray2d& C);
	
	void H_flux_1d_SD2(const doublearray1d& u_w, const doublearray1d& u_e, const double& gamma, const double& B1, doublearray1d& H);
	
	void end_of_step_1d_SD2(const doublearray2d& un, const double& dt, const double& t, double& dt_out, double& t_out, long& n, const double& gamma, const double& sum_t, const double& dt_cpu, const double& B1);
	
	double sign(const double& x);
	
	double min(const double& x, const double& y);
	
	double minmod(const double& x, const double& y);
	
	double minmod3(const double& x, const double& y, const double& z);
	
	void run_info_1d(double& dt, double& sum_t, long& J, double& cfl);

////////////////////////////////////////////////////////////////////////////////
// AUXILIARY FUNCTIONS -- User defined, see examples
////////////////////////////////////////////////////////////////////////////////

	void boundary_conditions(doublearray2d& u);

	void flux_x(const doublearray1d& v, doublearray1d& y, const double& gamma, const double& B1);

	void initial_conditions(doublearray2d& un, const double& x_init, const double& x_final, const double& gamma, doublearray1d& x);

	void spectral_radius(const doublearray1d& u, const double& gamma, const double& B1, double& rx);

	void writeout(const doublearray2d& un, const double& t, const double& gamma, const double& B1, const long& n);
}

#endif
