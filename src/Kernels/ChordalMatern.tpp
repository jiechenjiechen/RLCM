#ifndef _CHORDAL_MATERN_TPP_
#define _CHORDAL_MATERN_TPP_


//--------------------------------------------------------------------------
template<class Point>
double ChordalMatern::
Eval(const Point &x, const Point &y, double lambda) const {
  if (x.GetD() != 2 || y.GetD() != 2) {
    printf("ChordalMatern::Eval. Error: Point dimensions must be 2. Return NAN.\n");
    return NAN;
  }
  double r1 = x.Dist1(y);
  double r, ret;
  if (r1 == 0) {
    ret = s + lambda;
  }
  else {
    double *mx = x.GetPointer();
    double *my = y.GetPointer();
    double latx = mx[0], lonx = mx[1];
    double laty = my[0], lony = my[1];
    r = 2.0 * sqrt( Square(sin((latx-laty)/2.0)) +
                    cos(latx) * cos(laty) * Square(sin((lonx-lony)/2.0)) );
    r = sqrt2nu * r / ell;
    ret = s * pow(r,nu) * besselk(nu,r) / k0;
  }
  return ret;
}


#endif
