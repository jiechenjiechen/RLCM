#ifndef _ANISOTROPIC_GAUSSIAN_TPP_
#define _ANISOTROPIC_GAUSSIAN_TPP_


//--------------------------------------------------------------------------
template<class Point>
double AnisotropicGaussian::
Eval(const Point &x, const Point &y, double lambda) const {
  if (x.GetD() != y.GetD()) {
    printf("AnisotropicGaussian::Eval. Error: Point dimensions do not match. Return NAN.\n");
    return NAN;
  }
  if (x.GetD() != d) {
    printf("AnisotropicGaussian::Eval. Error: Points are not %ld dimensional. Return NAN.\n", (long)d);
    return NAN;
  }
  double r2 = x.Dist2(y, sigma);
  double ret = s * exp(-r2/2.0);
  if (r2 == 0.0) {
    ret += lambda;
  }
  return ret;
}


#endif
