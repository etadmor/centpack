#include "centpack_1d_FD2.h"

using namespace std;

void CENTPACK::run_info_1d(double& dt, double& sum_t, long& J, double& cfl)
{
	ofstream OutFile;
	
	OutFile.open("run_info.txt");
	OutFile<<"dt = "<<dt<<endl<<"total time = "<<sum_t<<endl<<"grid size = "<<J<<endl<<"CFL = "<<cfl<<endl;
	OutFile.close();
}
