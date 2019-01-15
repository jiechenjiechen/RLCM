#ifndef _ISOTROPIC_LAPLACE_TPP_
#define _ISOTROPIC_LAPLACE_TPP_


//--------------------------------------------------------------------------
template<class Point>
double IsotropicLaplace::
Eval(const Point &x, const Point &y, double lambda) const {
  if (x.GetD() != y.GetD()) {
    printf("IsotropicLaplace::Eval. Error: Point dimensions do not match. Return NAN.\n");
    return NAN;
  }
  double r2 = x.Dist2(y);
  double ret = s * exp(-sqrt(r2)/sigma);
  if (r2 == 0.0) {
    ret += lambda;
  }
  return ret;
}


#endif
