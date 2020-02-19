#include "centpack_2d_FD2.h"
#include "centpack_2d_SD2.h"
#include "centpack_2d_SD3.h"

using namespace std;

void CENTPACK::run_info_2d(double& dt, double& sum_t, long& J, long& K, double& cfl)
{
	ofstream OutFile;
	
	OutFile.open("run_info.txt");
	OutFile<<"dt = "<<dt<<endl<<"total time = "<<sum_t<<endl<<"grid size = "<<J<<" x "<<K<<endl<<"CFL = "<<cfl<<endl;
	OutFile.close();
}
