#ifndef _PCG_TPP_
#define _PCG_TPP_


//--------------------------------------------------------------------------
template<class MatrixA, class MatrixM>
void PCG::
Solve(MatrixA &A,        // Martix A
      MatState mStateA,  // Whether A is factorized
      const DVector &b,  // Right-hand side b
      const DVector &x0, // Initial guess x0
      MatrixM &M,        // Preconditioner M approx inv(A)
      MatState mStateM,  // Whether M is factorized
      INTEGER MaxIt,     // Maximum # of iterations
      double RTol,       // Relative residual tolerance
      bool TrueRes       // If true, use the true residual
                         // vector rathar than r to monitor
                         // convergence
      ) {

  NormB = b.Norm2();
  double Tol = RTol * NormB;
  Delete_1D_Array<double>(&mRes);
  New_1D_Array<double, INTEGER>(&mRes, MaxIt);

  mIter = 0;

  // r = b - Ax
  x = x0;
  DVector Ax;
  DVector r;
  if (x0.Norm1() == 0) {
    r = b;
  }
  else {
    A.MatVec(x, Ax, NORMAL, mStateA);
    b.Subtract(Ax, r);
  }
  mRes[mIter++] = r.Norm2();

  // true r
  DVector true_r;
  if (TrueRes) {
    true_r = r;
  }

  if (mRes[0] < Tol) {
    return;
  }

  // z = M(r)
  DVector z;
  M.MatVec(r, z, NORMAL, mStateM);

  // p = z
  DVector p = z;

  // rz = r'*z
  double rz = r.InProd(z);

  while (mIter < MaxIt) {

    // Ap = A(p)
    DVector Ap;
    A.MatVec(p, Ap, NORMAL, mStateA);

    // alpha = rz / (Ap'*p)
    double App = Ap.InProd(p);
    double alpha = rz / App;

    // x = x + alpha*p
    DVector v;
    p.Multiply(alpha, v);
    x.Add(v);

    // r = r - alpha*Ap
    Ap.Multiply(alpha, v);
    r.Subtract(v);

    // true r
    if (TrueRes) {
      A.MatVec(x, Ax, NORMAL, mStateA);
      b.Subtract(Ax, true_r);
      mRes[mIter] = true_r.Norm2();
    }
    else {
      mRes[mIter] = r.Norm2();
    }

    if (mRes[mIter] < Tol) {
      mIter++;
      break;
    }

    // z = M(r)
    M.MatVec(r, z, NORMAL, mStateM);

    // rz_new = r'*z
    double rz_new = r.InProd(z);

    // beta = rz_new / rz
    double beta = rz_new / rz;

    // rz = rz_new
    rz = rz_new;

    // p = z + beta*p
    p.Multiply(beta, v);
    z.Add(v, p);

    mIter++;

  }

}


#endif
