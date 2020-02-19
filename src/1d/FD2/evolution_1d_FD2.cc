#include "centpack_1d_FD2.h"

using namespace std;

void CENTPACK::evolution_1d_FD2(doublearray2d& un, const doublearray2d& uj_half, const double& lambda, const double& dx, const double& gamma, const double& alpha, const double& B1, const bool& odd)
{

	long j, l, J, L;
	
	J=un.getIndex1Size() - 4;
	L=un.getIndex2Size();
	
	doublearray2d un_half(J+4,L);
	
	predictor_1d_FD2(un, un_half, lambda, alpha, gamma, B1);
	corrector_1d_FD2(un, uj_half, un_half, lambda, gamma, B1, odd);
	
	boundary_conditions(un, odd);
}
