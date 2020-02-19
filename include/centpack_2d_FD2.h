#ifndef CENTPACK_2d_FD2_H
#define CENTPACK_2d_FD2_H

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
	
	int centpack_main_2d_FD2();

	void disclaimer();
	
	void read_in_2d(double& x_init, double& x_final, double& y_init, double& y_final, long& J, long& K, long& L, double& gamma, double& t_final, double& dt_out, double& cfl, double& alpha);
	
	void time_step_2d(const doublearray3d& un, const double& dx, const double& dy, const double& cfl, double& dt, double& t, double& dt_out, double& t_out, double& lambda, double& mu, const double& gamma);
	
	////////////////////////////////////////////////////////////////////////////
	// RECONSTRUCTION
	////////////////////////////////////////////////////////////////////////////
	
	void reconstruction_2d_FD2(const doublearray3d& un, doublearray3d& unhalf, const double& alpha, const bool& odd);

	// LIMITERS	-- minmod functions
	
	double sign(const double& x);

	double min(const double& x, const double& y);

	double minmod(const double& x, const double& y);

	double minmod3(const double& x, const double& y, const double& z);

	// end of reconstruction
	
	////////////////////////////////////////////////////////////////////////////
	// EVOLUTION
	////////////////////////////////////////////////////////////////////////////

	void evolution_2d_FD2(doublearray3d& un, const doublearray3d& unhalf, const double& lambda, const double& mu, const double& dx, const double& dy, const double& gamma, const double& alpha, const bool& odd);

	void predictor_2d_FD2(const doublearray3d& un, doublearray3d& u_star, const double& lambda, const double& mu, const double& dx, const double& dy, const double& gamma, const double& alpha);
	
	void corrector_2d_FD2(doublearray3d& un, const doublearray3d& unhalf, const doublearray3d& u_star, const double& lambda, const double& mu, const double& dx, const double& dy, const double& gamma, const double& alpha, const bool& odd);
	
	//end of evolution routines
	
	void end_of_step_2d_FD2(const doublearray3d& un, const double& dt, const double& t, double& dt_out, double& t_out, long& n, const double& gamma, const double& sum_t, const double& dt_cpu, bool& odd);

	void run_info_2d(double& dt, double& sum_t, long& J, long& K, double& cfl);
	
////////////////////////////////////////////////////////////////////////////////
// AUXILIARY FUNCTIONS -- User defined, see examples
////////////////////////////////////////////////////////////////////////////////

	void boundary_conditions(doublearray3d& u, const bool& odd);

	void flux_x(const doublearray1d& u, doublearray1d& f, const double& gamma);

	void flux_y(const doublearray1d& u, doublearray1d& g, const double& gamma);
	
	void initial_conditions(doublearray3d& un, const double& x_init, const double& x_final, const double& y_init, const double& y_final, double& dx, double& dy, const double& gamma, doublearray1d& x, doublearray1d& y);

	void spectral_radii(const doublearray1d& u, const double& gamma, double& rx, double& ry);

	void writeout(const doublearray3d& un, const double& t, const double& gamma, const long& n);
}

#endif
