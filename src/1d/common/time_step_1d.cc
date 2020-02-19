#include "centpack_1d_FD2.h"

using namespace std;

void CENTPACK::time_step_1d(const doublearray2d& un, const double& dx, const double& cfl, double& dt, double& t, double& dt_out, double& t_out, double& lambda, const double& gamma, const double& B1)
{
	long J=un.getIndex1Size() - 4;
	long L=un.getIndex2Size();
	long j, l;
	double rx, r_maxx;
	doublearray1d u_vector(L);
	
	r_maxx=0.0;
	
	for (j=2; j<J+2; j++)
	{
		for (l=0; l<L; l++)
			u_vector(l)=un(j,l);
		
		spectral_radius(u_vector, gamma, B1, rx);
		
		if (rx>r_maxx)
			r_maxx=rx;
	}
	
	dt=dx*cfl/r_maxx;
	dt = min(dt, dt_out - t_out);
	
	lambda=dt/dx;
	
	t=t+dt;
	t_out=t_out+dt;
}
