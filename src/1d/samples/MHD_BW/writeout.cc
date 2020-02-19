#include "centpack_1d_FD2.h"
#include "centpack_1d_SD2.h"
#include "centpack_1d_SD3.h"

using namespace std;
using std::string;

void CENTPACK::writeout(const doublearray2d& un, const double& t, const double& gamma, const double& B1, const long& n)
{
	
  cout.setf(ios::scientific, ios::floatfield);

  long j, J;
  char rho_file[20];
  char u1_file[20];
  char u2_file[20];
  char u3_file[20];
  char b2_file[20];
  char b3_file[20];
  char p_file[20];
  char t_file[20];
  
  J=un.getIndex1Size() - 4;
  
  doublearray1d rho(J);
  doublearray1d u1(J);
  doublearray1d u2(J);
  doublearray1d u3(J);
  doublearray1d b2(J);
  doublearray1d b3(J);
  doublearray1d p(J);
  
  for (j=0; j<J; j++)
  {
	  rho(j) = un(j+2,0);
	  u1(j) = un(j+2,1)/un(j+2,0);
	  u2(j) = un(j+2,2)/un(j+2,0);
	  u3(j) = un(j+2,3)/un(j+2,0);
	  b2(j) = un(j+2,4);
	  b3(j) = un(j+2,5);
	  p(j) = un(j+2,6) - .5*((pow(un(j+2,1),2.0) + pow(un(j+2,2),2.0) + pow(un(j+2,3),2.0))/un(j+2,0)) - .5*(pow(B1,2.0) + pow(un(j+2,4),2.0) + pow(un(j+2,5),2.0));
  }

  sprintf(rho_file, "rho_files/rho_%d", n);
  sprintf(u1_file, "u1_files/u1_%d", n);
  sprintf(u2_file, "u2_files/u2_%d", n);
  sprintf(u3_file, "u3_files/u3_%d", n);
  sprintf(b2_file, "b2_files/b2_%d", n);
  sprintf(b3_file, "b3_files/b3_%d", n);
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
  
  OutFile.open(u3_file, ios::out);
  OutFile<<u3;
  OutFile.close();
  
  OutFile.open(b2_file, ios::out);
  OutFile<<b2;
  OutFile.close();
  
  OutFile.open(b3_file, ios::out);
  OutFile<<b3;
  OutFile.close();

  OutFile.open(p_file, ios::out);
  OutFile<<p;
  OutFile.close();
  
  OutFile.open(t_file, ios::out);
  OutFile<<t;
  OutFile.close();
}
