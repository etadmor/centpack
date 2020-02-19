#include "centpack_1d_SD3.h"

using namespace std;

void CENTPACK::reconstruction_1d_SD3(doublearray2d& un, doublearray2d& u_E, doublearray2d& u_W,  const doublearray1d& ISl, const doublearray1d& ISc, const doublearray1d& ISr)
{

	long J=un.getIndex1Size() - 4;
	long L=un.getIndex2Size();
	long j, l;
	
	doublearray1d pl(2), pc(3), pr(2);
	
	double wl, wc, wr; 
	double alpl, alpc, alpr, alp_sum;
	double cl=.25, cc=.5, cr=.25;
	double eps=0.000001;
		
	for (l=0; l<L; l++)
	{		
		for (j=2; j<J+2; j++)
		{
			pl(0)=un(j,l);
			pl(1)=un(j,l)-un(j-1,l);
	
			pr(0)=pl(0);
			pr(1)=un(j+1,l)-un(j,l);
	
			pc(0)=pl(0) - (1.0/12.0)*(pr(1)-pl(1));
			pc(1)=.5*(pl(1)+pr(1));
			pc(2)=pr(1)-pl(1);
				
			alpl=cl/((eps+ISl(j))*(eps+ISl(j)));
			alpc=cc/((eps+ISc(j))*(eps+ISc(j)));
			alpr=cr/((eps+ISr(j))*(eps+ISr(j)));
				
			alp_sum=alpl+alpc+alpr;
				
			wl=alpl/alp_sum;
			wc=alpc/alp_sum;
			wr=alpr/alp_sum;
	
			u_E(j,l)=wl*(pl(0)+.5*pl(1)) + wc*(pc(0)+.5*pc(1)+.25*pc(2)) +wr*(pr(0)+.5*pr(1));
			u_W(j,l)=wl*(pl(0)-.5*pl(1)) + wc*(pc(0)-.5*pc(1)+.25*pc(2)) +wr*(pr(0)-.5*pr(1));
		}
	}
}
