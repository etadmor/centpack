#include "centpack_1d_FD2.h"
#include "centpack_1d_SD2.h"
#include "centpack_1d_SD3.h"

using namespace std;

void CENTPACK::boundary_conditions(doublearray2d& u)
{

  long J=u.getIndex1Size() - 4;
  long L=u.getIndex2Size();
  long j, l;
  
  for (l=0; l<L; l++)
  {
	  u(0,l) = 1.0;
	  u(1,l) = 1.0;
	  u(J+2,l) = 0.0;
	  u(J+3,l) = 0.0;
  }
}
