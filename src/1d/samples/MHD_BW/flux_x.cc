#include "centpack_1d_FD2.h"
#include "centpack_1d_SD2.h"
#include "centpack_1d_SD3.h"

using namespace std;

void CENTPACK::flux_x(const doublearray1d& v, doublearray1d& y, const double& gamma, const double& B1)
{
  double p, p_star;

  p=v(6) - .5*((pow(v(1),2.0) + pow(v(2),2.0) + pow(v(3),2.0))/v(0)) - .5*(pow(B1,2.0) + pow(v(4),2.0) + pow(v(5),2.0));
  p_star=p + .5*(pow(B1,2.0) + pow(v(4),2.0) + pow(v(5),2.0));

  y(0)=v(1);
  y(1)=(pow(v(1),2.0)/v(0)) + p_star;
  y(2)=(v(1)*v(2)/v(0)) - B1*v(4);
  y(3)=(v(1)*v(3)/v(0)) - B1*v(5);
  y(4)=(v(4)*v(1)/v(0)) - B1*v(2)/v(0);
  y(5)=(v(5)*v(1)/v(0)) - B1*v(3)/v(0);
  y(6)=(v(6) + p_star)*(v(1)/v(0)) - B1*(B1*(v(1)/v(0)) + (v(2)/v(0))*v(4) + (v(3)/v(0))*v(5));
}




















