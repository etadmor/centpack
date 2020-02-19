#include "centpack_1d_FD2.h"
#include "centpack_1d_SD2.h"
#include "centpack_1d_SD3.h"

using namespace std;

void CENTPACK::boundary_conditions(doublearray2d& u)
{

  long Nu=u.getIndex1Size() - 4;
  long L=u.getIndex2Size();
  long nu, l;
  
  for (l=0; l<L; l++)
  {
	  u(0,l) = u(Nu,l);
	  u(1,l) = u(Nu+1,l);
	  u(Nu+2,l) = u(2,l);
	  u(Nu+3,l) = u(3,l);
  }
}
