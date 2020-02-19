////////////////////////////////////////////////////////////////////////////////
//
// write_out.cc -- an auxiliary function of
//
// Central Hyperbolic Solver (C) Jorge Balbas and Eitan Tadmor, Nov. 2004
//
// Function called by solver.cc (main)
//
// This function outputs files containing the conserved solution variables (or 
// corresponding primitive variables) for the ideal MHD model.
// 
// Input (passed by reference):
//
// (1) un  --  a doublearray3d type variable holding the cell averages of u
//             at time t over the discretized solution domain
//
// (2) gamma -- a double type variable holding the ratio of specific heats
//
// (3) t -- time of solution being output
//
// Output (returned by reference):
//
// (1) Nine text files holding the solution variables at time t (n keeps track 
//     of how many output calls have been executed over the time of the 
//     simulation
//
//     rho_n -- density
//
//     u1_n -- velocity in x-direction
//
//     u2_n -- velocity in y-direction
//
//     p_n -- pressure
//
//     t_n -- time
//
////////////////////////////////////////////////////////////////////////////////

#include "centpack_2d_FD2.h"
#include "centpack_2d_SD2.h"
#include "centpack_2d_SD3.h"

using namespace std;

void CENTPACK::writeout(const doublearray3d& un, const double& t, const double& gamma, const long& n)
{

  cout.setf(ios::scientific, ios::floatfield);

  long j, k, J, K;
  char rho_file[20];
  char u1_file[20];
  char u2_file[20];
  char p_file[20];
  char t_file[20];
  
  J=un.getIndex1Size() - 4;
  K=un.getIndex2Size() - 4;
  
  doublearray2d rho(J,K);
  doublearray2d u1(J,K);
  doublearray2d u2(J,K);
  doublearray2d p(J,K);
  
  for (j=0; j<J; j++)
  {
	  for (k=0; k<K; k++)
	  {
		  rho(j,k)=un(j+2,k+2,0);
		  u1(j,k)=un(j+2,k+2,1)/un(j+2,k+2,0);
		  u2(j,k)=un(j+2,k+2,2)/un(j+2,k+2,0);
		  p(j,k)=(gamma-1.0)*(un(j+2,k+2,3)-.5*(pow(un(j+2,k+2,1),2) + pow(un(j+2,k+2,2),2))/un(j+2,k+2,0));
	  }
  }

  sprintf(rho_file, "rho_files/rho_%d", n);
  sprintf(u1_file, "u1_files/u1_%d", n);
  sprintf(u2_file, "u2_files/u2_%d", n);
  sprintf(p_file, "p_files/p_%d", n);
  sprintf(t_file, "t_files/t_%d", n);
  
  ofstream OutFile;
  OutFile.open(rho_file, ios::out);
  OutFile<<rho;
  OutFile.close();
  
  OutFile.open(u1_file, ios::out);
  OutFile<<u1;
  OutFile.close();
  
  OutFile.open(u2_file, ios::out);
  OutFile<<u2;
  OutFile.close();

  OutFile.open(p_file, ios::out);
  OutFile<<p;
  OutFile.close();
  
  OutFile.open(t_file, ios::out);
  OutFile<<t;
  OutFile.close();
}
