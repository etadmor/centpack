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
  char t_file[20];
  
  J=un.getIndex1Size() - 4;
  
  doublearray1d rho(J);
  
  for (j=0; j<J; j++)
	  rho(j) = un(j+2,0);

  sprintf(rho_file, "rho_files/rho_%d", n);
  sprintf(t_file, "t_files/t_%d", n);
  
  ofstream OutFile;
  OutFile.open(rho_file, ios::out);
  OutFile<<rho;
  OutFile.close();
  
  OutFile.open(t_file, ios::out);
  OutFile<<t;
  OutFile.close();
}
