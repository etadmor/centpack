#include "centpack_1d_FD2.h"

using namespace std;

void CENTPACK::end_of_step_1d_FD2(const doublearray2d& un, const double& dt, const double& t, double& dt_out, double& t_out, long& n, const double& gamma, const double& sum_t, const double& dt_cpu, const double& B1, const bool& odd)
{
	if (t_out == dt_out)
	{
		writeout(un, t, gamma, B1, n);
		n = n++;
		
		t_out = 0.0;
		
		cout<<"output written at t = "<<t<<endl;
	}
	
	cout<<"run="<<t<<",dt="<<dt<<",cpu_t="<<dt_cpu<<",t="<<sum_t<<",odd = "<<odd<<endl;
}
