#include "centpack_1d_FD2.h"
#include "centpack_1d_SD2.h"
#include "centpack_1d_SD3.h"

using namespace std;

void CENTPACK::initial_conditions(doublearray2d& un, const double& x_init, const double& x_final, const double& gamma, doublearray1d& x)
{

  long J=un.getIndex1Size() - 4;
  long L=un.getIndex2Size();
  long j, l;
  long j1 = J/2 + 2;
  double dx;
 
  doublearray1d u_left(L), u_right(L);
  
  dx=(x_final - x_init)/J;
  
  x(0)=x_init - 2.0*dx;
  
  for (j=1; j<J+4; j++)
	  x(j)=x(j-1)+dx;
  
  u_left(0)=1.0;
  u_left(1)=0.0;
  u_left(2)=0.0;
  u_left(3)=0.0;
  u_left(4)=1.0;
  u_left(5)=0.0;
  u_left(6)=1.0+(25.0/32.0);

  u_right(0)=.125; 
  u_right(1)=0.0;
  u_right(2)=0.0;
  u_right(3)=0.0;
  u_right(4)=-1.0;
  u_right(5)=0.0;
  u_right(6)=.1+(25.0/32.0);
  
  for(l=0; l<L; l++)
  {
	  for (j=0; j<j1; j++)
		  un(j,l)=u_left(l);
	  for (j=j1; j<J+4; j++)
		  un(j,l)=u_right(l);
  }
}
