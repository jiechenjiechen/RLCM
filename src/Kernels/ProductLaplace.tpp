#ifndef _PRODUCT_LAPLACE_TPP_
#define _PRODUCT_LAPLACE_TPP_


//--------------------------------------------------------------------------
template<class Point>
double ProductLaplace::
Eval(const Point &x, const Point &y, double lambda) const {
  if (x.GetD() != y.GetD()) {
    printf("ProductLaplace::Eval. Error: Point dimensions do not match. Return NAN.\n");
    return NAN;
  }
  double r1 = x.Dist1(y);
  double ret = s * exp(-r1/sigma);
  if (r1 == 0.0) {
    ret += lambda;
  }
  return ret;
}


#endif
