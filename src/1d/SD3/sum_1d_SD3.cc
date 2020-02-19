#include "centpack_1d_SD3.h"

using namespace std;

double CENTPACK::sum_1d_SD3(const doublearray1d& u)
{
	long L=u.getIndex1Size();
	long l;
	double s=0.0;
	
	for (l=0; l<L; l++)
		s += u(l);
	
	return s;
}
