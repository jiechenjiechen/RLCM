#ifndef _INV_MULTIQUADRIC_TPP_
#define _INV_MULTIQUADRIC_TPP_


//--------------------------------------------------------------------------
template<class Point>
double InvMultiquadric::
Eval(const Point &x, const Point &y, double lambda) const {
  if (x.GetD() != y.GetD()) {
    printf("InvMultiquadric::Eval. Error: Point dimensions do not match. Return NAN.\n");
    return NAN;
  }
  double r2 = x.Dist2(y);
  double ret = s / sqrt( r2/Square(sigma) + 1.0 );
  if (r2 == 0.0) {
    ret += lambda;
  }
  return ret;
}


#endif
