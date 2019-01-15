#ifndef _POLYNOMIAL_TPP_
#define _POLYNOMIAL_TPP_


//--------------------------------------------------------------------------
template<class Point>
double Polynomial::
Eval(const Point &x, const Point &y, double lambda) const {
  if (x.GetD() != y.GetD()) {
    printf("Polynomial::Eval. Error: Point dimensions do not match. Return NAN.\n");
    return NAN;
  }
  double r1 = x.Dist1(y);
  double ret = s * pow(a * x.InProd(y) + c, deg);
  if (r1 == 0.0) {
    ret += lambda;
  }
  return ret;
}


#endif
