////////////////////////////////////////////////////////////////////////////////
//
// CentPack -- A generic numerical solver for hyperbolic conservation
//             laws and related time dependent problems
//
// Copyright (C) 2006 Jorge Balbas and Eitan Tadmor
// 
// initial_conditions.cc -- an auxiliary function of
//
// Function called by cenatpack_main, problem specific.
//
// This function initializes the solution array with the initial conditions of 
// the problem to be solved
// 
// Input (passed by reference):
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
// (5) dx -- space scale in x-direction
//
// (6) dz -- space scale in z-direction
//
// (7) gamma - a double type variable holding the ratio of specific heats
//
// Output (returned by reference):
//
// (1) un -- a doublearray3d type variable holding the conserved
//           quantities density, momentum (3 components), magnetic field
//           (3 components), and total energy over the entire soltion domain
//           at t=0
//
// (2) x -- a doublearray1d type variable holding the discretization over the
//          interval [x_init, x_final]
//
// (3) z -- a doublearray1d type variable holding the discretization over the
//          interval [z_init, z_final]
//
////////////////////////////////////////////////////////////////////////////////

#include "centpack_2d_FD2.h"
#include "centpack_2d_SD2.h"
#include "centpack_2d_SD3.h"

using namespace std;

void CENTPACK::initial_conditions(doublearray3d& un, const double& x_init, const double& x_final, const double& y_init, const double& y_final, double& dx, double& dy, const double& gamma, doublearray1d& x, doublearray1d& y)
{

  long J=un.getIndex1Size() - 4;
  long K=un.getIndex2Size() - 4;
  long L=un.getIndex3Size();
  long j, k, l;
  double I1, I2, I3, I4;
  
  dx=(x_final - x_init)/J;
  dy=(y_final - y_init)/K;
  
  x(0)=x_init - 2.0*dx;
  y(0)=y_init - 2.0*dy;
  
  for (j=1; j<J+4; j++)
	  x(j)=x(j-1)+dx;
  
  for (k=1; k<K+4; k++)
	  y(k)=y(k-1)+dy;

  for (k=0; k<K+4; k++)
  {
	  for (j=0; j<J+4; j++)
	  {
		  un(j,k,0)=pow(gamma,2.0);
		  un(j,k,1)=un(j,k,0)*(1.0/dy)*(cos(y(k)+.5*dy)-cos(y(k)-.5*dy));
		  un(j,k,2)=0.0;
		  un(j,k,3)=-un(j,k,0)*(1.0/dx)*(cos(x(j)+.5*dx)-cos(x(j)-.5*dx));
		  un(j,k,4)=(1.0/dy)*(cos(y(k)+.5*dy)-cos(y(k)-.5*dy));
		  un(j,k,5)=0.0;
		  un(j,k,6)=(-1.0/(2*dx))*(cos(2*(x(j)+.5*dx))-cos(2*(x(j)-.5*dx)));
		  
		  I1=-(0.125/dy)*(un(j,k,0)+1.0)*(sin(2*(y(k)+.5*dy)) - sin(2*(y(k)-.5*dy)));
		  I2=-(0.125/dx)*un(j,k,0)*(sin(2*(x(j)+.5*dx)) - sin(2*(x(j)-.5*dx)));
		  I3=-(0.0625/dx)*(sin(4*(x(j)+.5*dx)) - sin(4*(x(j)-.5*dx)));
		  un(j,k,7)=3.0 + .5*un(j,k,0) + I1 + I2 + I3;
	}
  }
}
