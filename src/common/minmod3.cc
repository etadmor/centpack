#include "centpack_1d_FD2.h"

using namespace std;

double CENTPACK::sign(const double& x)
{
  if (x<0)
    return -1.0;
  else
    return 1.0;
}

double CENTPACK::min(const double& x, const double& y)
{
	return x < y ? x : y ;
}

double CENTPACK::minmod(const double& x, const double& y)
{
	return .5*(sign(x) + sign(y))*min(fabs(x),fabs(y));
}

double CENTPACK::minmod3(const double& x, const double& y, const double& z)
{    
	return minmod(x,minmod(y,z));
}
