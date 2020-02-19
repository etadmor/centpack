#include "centpack_1d_SD3.h"

using namespace std;

int CENTPACK::centpack_main_1d_SD3()
{
	disclaimer();
	
	double x_init, x_final, t, t_init, t_final, t_out;
	double B1;
	double dx, dt, cfl, dt_out, alpha;
	double lambda;
	double gamma;
	double sum_t=0.0;
	double dt_cpu=0.0;
	double t_start=clock();
	long J, L;
	long j, l;
	long n=0;
	
	read_in_1d(x_init, x_final, J, L, gamma, t_final, dt_out, cfl, alpha, B1);
	
	doublearray1d x(J+4);
	
	cout.setf(ios::scientific, ios::floatfield);
	
	t=0.0;
	t_init=0.0;
	dx=(x_final - x_init)/J;
	
	doublearray2d un(J+4,L);
	
	initial_conditions(un, x_init, x_final, gamma, x);
	writeout(un, t, gamma, B1, n);
	
	t_out=0.0;
	n++;
	
	do
	{
		time_step_1d(un, dx, cfl, dt, t, dt_out, t_out, lambda, gamma, B1);
		evolution_1d_SD3(un, lambda, dx, gamma, B1);
		
		dt_cpu=(clock()-t_start)/CLOCKS_PER_SEC;
		sum_t=sum_t+dt_cpu;
		
		end_of_step_1d_SD3(un, dt, t, dt_out, t_out, n, gamma, sum_t, dt_cpu, B1);
		
		t_start=clock();
			
	}while(t<t_final);
	
	writeout(un, t, gamma, B1, n);
	
	run_info_1d(dt, sum_t, J, cfl);
	
	return 0;

}
//////////////////////////////////////////////////////////////////////////////
