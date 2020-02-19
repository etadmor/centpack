////////////////////////////////////////////////////////////////////////////
//
// CentPack -- A generic numerical solver for hyperbolic conservation
//             laws and related time dependent problems
//
// Copyright (C) 2006 Jorge Balbas and Eitan Tadmor
//
// read_in.cc -- reads input
//
// Called by solver (main)
// 
// This function reads input parameters for the simulation provided by the user // in a separate file named "input"
//
// Input (passed by reference):
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
// This function does not produce any direct output, it only initializes the
// input variables described above with the values provided by the user
//
////////////////////////////////////////////////////////////////////////////////

#include "centpack_2d_FD2.h"
#include "centpack_2d_SD2.h"
#include "centpack_2d_SD3.h"

using namespace std;

void CENTPACK::read_in_2d(double& x_init, double& x_final, double& y_init, double& y_final, long& J, long& K, long& L, double& gamma, double& t_final, double& dt_out, double& cfl, double& alpha)
{	
	ifstream InFile;
	
	InFile.open("input", ios::in);
	
	if (!InFile)
	{
		cerr<<"Error in opening input file";
		exit(1);
	}
	
	InFile>>x_init;
	InFile>>x_final;
	InFile>>y_init;
	InFile>>y_final;
	InFile>>J;
	InFile>>K;
	InFile>>L;
	InFile>>gamma;
	InFile>>t_final;
	InFile>>dt_out;
	InFile>>cfl;
	InFile>>alpha;
}

