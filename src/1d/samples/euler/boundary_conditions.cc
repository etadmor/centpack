#include "centpack_1d_FD2.h"
#include "centpack_1d_SD2.h"
#include "centpack_1d_SD3.h"

using namespace std;

void CENTPACK::boundary_conditions(doublearray2d& u)
{

  long Nu = u.getIndex1Size() - 4;
  long L = u.getIndex2Size();
  long nu, l;
  
  double p_left;
  double p_right;
  double gamma = 1.4;
  
  doublearray1d u_left(L), u_right(L);
  
  p_left = 1.0;
  p_right = 0.1;
  
  u_left(0) = 1.0;
  u_left(1) = 0.0;
  u_left(2) = p_left/(gamma - 1.0);

  u_right(0) = 0.125; 
  u_right(1) = 0.0;
  u_right(2) = p_right/(gamma - 1.0);
  
  for (l = 0; l < L; l++)
  {
	  u(0,l) = u_left(l);
	  u(1,l) = u_left(l);
	  u(Nu+2,l) = u_right(l);
	  u(Nu+3,l) = u_right(l);
  }
}
