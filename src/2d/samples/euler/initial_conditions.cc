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

void CENTPACK::initial_conditions(doublearray3d& un, const double& x_init, const double& x_final, const double& y_init, const double& y_final,  double& dx, double& dy, const double& gamma, doublearray1d& x, doublearray1d& y)
{

  long J=un.getIndex1Size() - 4;
  long K=un.getIndex2Size() - 4;
  long L=un.getIndex3Size();
  long j, k, l;
  long j1 = J/2 + 2;
  long k1 = K/2 + 2;
  
  double p_one=1.5;
  double p_two=0.3;
  double p_three=0.029;
  double p_four=0.3;
  
  doublearray1d u_one(L), u_two(L), u_three(L), u_four(L);
  
  dx=(x_final - x_init)/J;
  dy=(y_final-y_init)/K;
  
  x(0)=x_init - 2.0*dx;
  y(0)=y_init - 2.0*dy;
  
  for (j=1; j<J+4; j++)
	  x(j)=x(j-1)+dx;
  
  for (k=1; k<K+4; k++)
	  y(k)=y(k-1)+dy;

  u_one(0)=1.5;
  u_one(1)=0.0;
  u_one(2)=0.0;
  u_one(3)=p_one/(gamma - 1.0) + .5*(pow(u_one(1),2.0) + pow(u_one(2),2.0))/u_one(0);
  
  u_two(0)=0.5323;
  u_two(1)=1.206*u_two(0);
  u_two(2)=0.0;
  u_two(3)=p_two/(gamma - 1.0) + .5*(pow(u_two(1),2.0) + pow(u_two(2),2.0))/u_two(0);
  
  u_three(0)=0.138;
  u_three(1)=1.206*u_three(0);
  u_three(2)=1.206*u_three(0);
  u_three(3)=p_three/(gamma - 1.0) + .5*(pow(u_three(1),2.0) + pow(u_three(2),2.0))/u_three(0);
  
  u_four(0)=0.5323;
  u_four(1)=0.0;
  u_four(2)=1.206*u_four(0);
  u_four(3)=p_four/(gamma - 1.0) + .5*(pow(u_four(1),2.0) + pow(u_four(2),2.0))/u_four(0);
  
  for(l=0; l<L; l++)
  {
	  for (k=k1; k<K+4; k++)
	  {
		  for (j=j1; j<J+4; j++)
			  un(j,k,l)=u_one(l);
		  for (j=0; j<j1; j++)
			  un(j,k,l)=u_two(l);
	  }
	  
	  for (k=0; k<k1; k++)
	  {
		  for (j=j1; j<J+4; j++)
			  un(j,k,l)=u_four(l);
		  for (j=0; j<j1; j++)
			  un(j,k,l)=u_three(l);
	  }
  }
}
