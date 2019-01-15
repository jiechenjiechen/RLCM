#ifndef _CHI2_TPP_
#define _CHI2_TPP_


//--------------------------------------------------------------------------
template<class Point>
double Chi2::
Eval(const Point &x, const Point &y, double lambda) const {
  if (x.GetD() != y.GetD()) {
    printf("Chi2::Eval. Error: Point dimensions do not match. Return NAN.\n");
    return NAN;
  }
  double r1 = x.Dist1(y);
  double ret = s * x.KernelFuncChi2(y);
  if (r1 == 0.0) {
    ret += lambda;
  }
  return ret;
}


#endif
