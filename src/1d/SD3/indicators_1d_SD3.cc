#include "centpack_1d_SD3.h"

using namespace std;

void CENTPACK::indicators_1d_SD3(doublearray2d& un, doublearray1d& ISl, doublearray1d& ISc, doublearray1d& ISr, const doublearray1d& u_norm)
{

	long J=un.getIndex1Size() - 4;
	long L=un.getIndex2Size();
	long j, l;
	
	doublearray1d pl(2), pc(3), pr(2);
	doublearray1d ISl_vector(L), ISc_vector(L), ISr_vector(L);
	double eps=0.000001;
	
	for (j=2; j<J+2; j++)
	{
		for (l=0; l<L; l++)
		{
			pl(0)=un(j,l);
			pl(1)=un(j,l)-un(j-1,l);
	
			pr(0)=pl(0);
			pr(1)=un(j+1,l)-un(j,l);
	
			pc(0)=pl(0) - (1.0/12.0)*(pr(1)-pl(1));
			pc(1)=.5*(pl(1)+pr(1));
			pc(2)=pr(1)-pl(1);
			
			ISl_vector(l)=pl(1)*pl(1)/(u_norm(l)+eps);
			ISc_vector(l)=(1.0/(u_norm(l)+eps))*((13.0/3.0)*pc(2)*pc(2) + pc(1)*pc(1));
			ISr_vector(l)=pr(1)*pr(1)/(u_norm(l)+eps);
		}
			
		ISl(j)=(1.0/L)*sum_1d_SD3(ISl_vector);
		ISc(j)=(1.0/L)*sum_1d_SD3(ISc_vector);
		ISr(j)=(1.0/L)*sum_1d_SD3(ISr_vector);
	}
}
