#include "centpack_1d_FD2.h"
#include "centpack_1d_SD2.h"
#include "centpack_1d_SD3.h"

using namespace std;
using std::string;

void CENTPACK::writeout(const doublearray2d& un, const double& t, const double& gamma, const double& B1, const long& n)
{
	
  cout.setf(ios::scientific, ios::floatfield);

  long j, J;
  char u_file[20];
  char t_file[20];
  
  J=un.getIndex1Size() - 4;
  
  doublearray1d u(J);
  
  for (j=0; j<J; j++)
	  u(j) = un(j+2,0);

  sprintf(u_file, "u_files/u_%d", n);
  sprintf(t_file, "t_files/t_%d", n);
  
  ofstream OutFile;
  OutFile.open(u_file, ios::out);
  OutFile<<u;
  OutFile.close();
  
  OutFile.open(t_file, ios::out);
  OutFile<<t;
  OutFile.close();
}
