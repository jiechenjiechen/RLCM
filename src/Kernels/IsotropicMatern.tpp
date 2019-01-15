#ifndef _ISOTROPIC_MATERN_TPP_
#define _ISOTROPIC_MATERN_TPP_


//--------------------------------------------------------------------------
template<class Point>
double IsotropicMatern::
Eval(const Point &x, const Point &y, double lambda) const {
  if (x.GetD() != y.GetD()) {
    printf("IsotropicMatern::Eval. Error: Point dimensions do not match. Return NAN.\n");
    return NAN;
  }
  double r2 = x.Dist2(y);
  double r, ret;
  if (r2 == 0) {
    ret = s + lambda;
  }
  else {
    r = sqrt2nu * sqrt(r2) / ell;
    ret = s * pow(r,nu) * besselk(nu,r) / k0;
  }
  return ret;
}


#endif
