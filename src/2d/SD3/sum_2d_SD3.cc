////////////////////////////////////////////////////////////////////////////
//
// CentPack -- A generic numerical solver for hyperbolic conservation
//             laws and related time dependent problems
//
// Copyright (C) 2006 Jorge Balbas and Eitan Tadmor
//
// sum.cc -- a core function of CentPack
//
// Function called by reconstruction.cc and reconstruccion_diag.cc
//
// This funcion calculates the sum of the entries of a doublearray1d type 
// variable
//
// Input (passed by reference):
//
// (1) u  --  a doublearray1d
//
// Output (returned by value):
//
// (1) s -- a double type variable holding sum of the entries of a doublearray1d
//
////////////////////////////////////////////////////////////////////////////////

#include "centpack_2d_SD3.h"

using namespace std;

double CENTPACK::sum_2d_SD3(const doublearray1d& u)
{
	long L=u.getIndex1Size();
	long l;
	double s=0.0;
	
	for (l=0; l<L; l++)
		s=s+u(l);
	
	return s;
}
