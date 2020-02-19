////////////////////////////////////////////////////////////////////////////////
//
// spctral_radii.cc -- an auxiliary function of
//
// CentPack (C) Jorge Balbas and Eitan Tadmor, 2006
//
// Function called by time_step_2d.cc, Hx_flux.cc and Hz_flux.cc (model 
// specific)
//
// This function calculates the spectral radii of the Jacobian matrices of
// f(u) and g(u) for ideal MHD equations
// 
// Input (passed by reference):
//
// (1) u -- a doublearray1d type variable holding the conserved
//          quantities density, momentum (3 components), magnetic field
//          (3 components), and total energy over the entire soltion domain
//          at a cell of the solution domain
//
// (2) gamma -- a double type variable holding the ratio of specific heats
//
// Output (returned by reference):
//
// (1) rx -- a double type variable holding the spectral radius of the 
//           Jacobian matrix of f(u)
//
// (2) rz -- a double type variable holding the spectral radius of the 
//           Jacobian matrix of g(u)
//
////////////////////////////////////////////////////////////////////////////////

#include "centpack_2d_FD2.h"
#include "centpack_2d_SD2.h"
#include "centpack_2d_SD3.h"

using namespace std;

void CENTPACK::spectral_radii(const doublearray1d& u, const double& gamma, double& rx, double& ry)
{
	double rho, vx, vy, p, H, c;
	
	rho=u(0);
	vx=u(1)/rho;
	vy=u(2)/rho;
	p=(gamma-1.0)*(u(3)-.5*rho*(pow(vx,2) + pow(vy,2)));
	c=sqrt(gamma*p/rho);
	rx=fabs(vx)+c;
	ry=fabs(vy)+c;
}
