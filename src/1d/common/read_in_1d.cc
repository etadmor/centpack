#include "centpack_1d_FD2.h"
#include "centpack_1d_SD2.h"

using namespace std;

void CENTPACK::read_in_1d(double& x_init, double& x_final, long& J, long& L, double& gamma, double& t_final, double& dt_out, double& cfl, double& alpha, double& B1)
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
	InFile>>J;
	InFile>>L;
	InFile>>gamma;
	InFile>>t_final;
	InFile>>dt_out;
	InFile>>cfl;
	InFile>>alpha;
	InFile>>B1;
}

