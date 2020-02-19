////////////////////////////////////////////////////////////////////////////////
//
// CentPack -- A generic numerical solver for hyperbolic conservation
//             laws and related time dependent problems
//
// Copyright (C) 2006 Jorge Balbas and Eitan Tadmor
// 
// boundary_conditions.cc
//
// Function called by evolution, problem specific.
//
// This function implements periodic boundary conditions at the four boundaries
// of a rectangular domain by filling the first and last two rows/columns of
// the solution array u according to the periodic conditions u(x,z)=u(x+p,z+p).
// 
// Input (passed by reference):
//
// (1)  u -- a doublearray3d type variable holding the conserved
//           quantities density, momentum (3 components), magnetic field
//           (3 components), and total energy over the entire soltion domain.
//
// Output (returned by reference):
//
// (1) u -- a doublearray3d type variable holding the conserved
//          quantities density, momentum (3 components), magnetic field
//          (3 components), and total energy over the entire soltion domain.
//
////////////////////////////////////////////////////////////////////////////////

#include"centpack_2d_FD2.h"
#include"centpack_2d_SD2.h"
#include"centpack_2d_SD3.h"

using namespace std;

void CENTPACK::boundary_conditions(doublearray3d& u)
{
  
  long J=u.getIndex1Size() - 4;
  long K=u.getIndex2Size() - 4;
  long L=u.getIndex3Size();
  long j, k, l;
  
  for (l=0; l<L; l++)
  {
	  for (j=2; j<J+2; j++)
	  {
		  u(j,0,l)=u(j,K,l);
		  u(j,1,l)=u(j,K+1,l);
		  u(j,K+2,l)=u(j,2,l);
		  u(j,K+3,l)=u(j,3,l);
	  }
	  
	  for (k=0; k<K+4; k++)
	  {
		  u(0,k,l)=u(J,k,l);
		  u(1,k,l)=u(J+1,k,l);
		  u(J+2,k,l)=u(2,k,l);
		  u(J+3,k,l)=u(3,k,l);
	  }
  }
}

void CENTPACK::boundary_conditions(doublearray3d& u, const bool& odd)
{
  
  long J=u.getIndex1Size() - 4;
  long K=u.getIndex2Size() - 4;
  long L=u.getIndex3Size();
  long j, k, l;
  
  if (odd)
  {
	  for (l=0; l<L; l++)
	  {
		  for (j=1; j<J+2; j++)
		  {
			  u(j,0,l)=u(j,K,l);
			  u(j,K+2,l)=u(j,2,l);
			  u(j,K+3,l)=u(j,3,l);
		  }
		  
		  for (k=0; k<K+4; k++)
		  {
			  u(0,k,l)=u(J,k,l);
			  u(J+2,k,l)=u(2,k,l);
			  u(J+3,k,l)=u(3,k,l);
		  }
	  }
  }
  
  else
  {
	  for (l=0; l<L; l++)
	  {
		  for (j=2; j<J+3; j++)
		  {
			  u(j,0,l)=u(j,K,l);
			  u(j,1,l)=u(j,K+1,l);
			  u(j,K+3,l)=u(j,3,l);
		  }
		  
		  for (k=0; k<K+4; k++)
		  {
			  u(0,k,l)=u(J,k,l);
			  u(1,k,l)=u(J+1,k,l);
			  u(J+3,k,l)=u(3,k,l);
		  }
	  }
  }
}
