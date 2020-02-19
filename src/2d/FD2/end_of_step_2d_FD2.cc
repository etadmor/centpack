#include "centpack_2d_FD2.h"

using namespace std;

void CENTPACK::end_of_step_2d_FD2(const doublearray3d& un, const double& dt, const double& t, double& dt_out, double& t_out, long& n, const double& gamma, const double& sum_t, const double& dt_cpu, bool& odd)
{
	if (t_out == dt_out)
	{
		writeout(un, t, gamma, n);
		n = n++;
		
		t_out = 0.0;
		
		cout<<"output written at t = "<<t<<endl;
	}
	
	cout<<"run="<<t<<",dt="<<dt<<",cpu_t="<<dt_cpu<<",t="<<sum_t<<",parity= "<<odd<<endl;
	
	odd = !(odd);
	
}
