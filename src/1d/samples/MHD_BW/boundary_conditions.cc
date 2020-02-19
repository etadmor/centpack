#include "centpack_1d_FD2.h"
#include "centpack_1d_SD2.h"
#include "centpack_1d_SD3.h"

using namespace std;

void CENTPACK::boundary_conditions(doublearray2d& u, const bool& odd)
{

  long J=u.getIndex1Size() - 4;
  long L=u.getIndex2Size();
  long j, l;
  
  doublearray1d u_left(L), u_right(L);
  
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
  
  if (odd)
  {
	  for (l=0; l<L; l++)
	  {
		  u(0,l) = u_left(l);
		  u(J+2,l) = u_right(l);
		  u(J+3,l) = u_right(l);
	  }
  }
  
  else
  {
	  for (l=0; l<L; l++)
	  {
		  u(0,l) = u_left(l);
		  u(1,l) = u_left(l);
		  u(J+3,l) = u_right(l);
	  }
  }
}

void CENTPACK::boundary_conditions(doublearray2d& u)
{

  long J=u.getIndex1Size() - 4;
  long L=u.getIndex2Size();
  long j, l;
  
  doublearray1d u_left(L), u_right(L);
  
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
  
  for (l=0; l<L; l++)
  {
	  u(0,l) = u_left(l);
	  u(1,l) = u_left(l);
	  u(J+2,l) = u_right(l);
	  u(J+3,l) = u_right(l);
  }
}
