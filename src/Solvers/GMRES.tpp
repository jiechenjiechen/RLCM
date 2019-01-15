#ifndef _GMRES_TPP_
#define _GMRES_TPP_


//--------------------------------------------------------------------------
template<class MatrixA, class MatrixM>
void GMRES::
Solve(MatrixA &A,        // Martix A
      MatState mStateA,  // Whether A is factorized
      const DVector &b,  // Right-hand side b
      const DVector &x0, // Initial guess x0
      MatrixM &M,        // Preconditioner M approx inv(A)
      MatState mStateM,  // Whether M is factorized
      INTEGER m,         // Restart cycle m
      INTEGER MaxIt,     // Maximum # of iterations
      double RTol        // Relative residual tolerance
      ) {
  SolveBase(A, mStateA, b, x0, M, mStateM, m, MaxIt, RTol, USUAL);
}


//--------------------------------------------------------------------------
template<class MatrixA, class MatrixM>
void GMRES::
RefineSolve(MatrixA &A,        // Martix A
            MatState mStateA,  // Whether A is factorized
            const DVector &b,  // Right-hand side b
            const DVector &x0, // Initial guess x0
            MatrixM &M,        // Preconditioner M approx inv(A)
            MatState mStateM   // Whether M is factorized
            ) {
  SolveBase(A, mStateA, b, x0, M, mStateM,
            MAX_REFINE_STEP, MAX_REFINE_STEP, EPS, REFINE);
}


//--------------------------------------------------------------------------
template<class MatrixA, class MatrixM>
void GMRES::
SolveBase(MatrixA &A,        // Martix A
          MatState mStateA,  // Whether A is factorized
          const DVector &b,  // Right-hand side b
          const DVector &x0, // Initial guess x0
          MatrixM &M,        // Preconditioner M approx inv(A)
          MatState mStateM,  // Whether M is factorized
          INTEGER m,         // Restart cycle m
          INTEGER MaxIt,     // Maximum # of iterations
          double RTol,       // Relative residual tolerance
          SolveType Type     // Normal solve or refinement solve?
          ) {

  // Norm of RHS, tolerance, residual history
  NormB = b.Norm2();
  double Tol = RTol * NormB;
  Delete_1D_Array<double>(&mRes);
  New_1D_Array<double, INTEGER>(&mRes, MaxIt);

  // Basis matrix V
  INTEGER N = b.GetN();
  DMatrix V(N, m+1);

  // The matrix H of size m*(m-1) is used to monitor convergence and
  // extract solution. For computational efficiency and ease of
  // operation, we abandom the use of a DMatrix object to store H, but
  // use the traditional double* array. Accompanied with the matrix H
  // are vector g and arrays c and s.
  double *H = NULL, *g = NULL, *c = NULL, *s = NULL;
  New_1D_Array<double, INTEGER>(&H, m*(m-1));
  New_1D_Array<double, INTEGER>(&g, m);
  New_1D_Array<double, INTEGER>(&c, m-1);
  New_1D_Array<double, INTEGER>(&s, m-1);

  // Iteration count
  mIter = 0;

  // r = b - Ax
  // beta = norm(r)
  // g[0] = beta
  x = x0;
  DVector Ax;
  DVector r;
  double beta;
  if (x0.Norm1() == 0) {
    r = b;
    beta = NormB;
  }
  else {
    A.MatVec(x, Ax, NORMAL, mStateA);
    b.Subtract(Ax, r);
    beta = r.Norm2();
  }
  g[0] = beta;
  mRes[mIter++] = beta;
  bool QuitIter = false;

  // Check convergence
  if (beta < Tol) {
    goto cleanup;
  }

  // Loop
  INTEGER j;
  while (mIter < MaxIt) {

    // V(:,0) = r / beta
    DVector v;
    r.Divide(beta, v);
    V.SetColumn(0, v);

    // For each restart cycle
    for (j = 1; j < m; j++) {

      // w = A(M(V(:,j-1)))
      DVector Mv, w;
      V.GetColumn(j-1, v);
      M.MatVec(v, Mv, NORMAL, mStateM);
      A.MatVec(Mv, w, NORMAL, mStateA);

      // Modified Gram-Schmidt. w = w - V(:,i)*(V(:,i)'*w)
      double hij;
      for (INTEGER i = 0; i < j; i++) {
        V.GetColumn(i, v);
        hij = v.InProd(w);
        H[i+(j-1)*m] = hij;
        DVector vv;
        v.Multiply(hij, vv);
        w.Subtract(vv);
      }

      // H(j,j-1) = norm(w)
      hij = w.Norm2();
      H[j+(j-1)*m] = hij;

      // V(:,j) = w / H(j,j-1)
      w.Divide(hij, v);
      V.SetColumn(j, v);

      // Solve (implicitly) the least squares problem
      // min(norm(beta*e_1-H*y)) and obtain residual
      for (INTEGER i = 0; i < j-1; i++) {
        drot_(&ONEi, H+i+(j-1)*m, &ONEi, H+(i+1)+(j-1)*m, &ONEi, c+i, s+i);
      }
      double da = H[(j-1)+(j-1)*m], db = H[j+(j-1)*m];
      drotg_(&da, &db, c+j-1, s+j-1);
      drot_(&ONEi, H+(j-1)+(j-1)*m, &m, H+j+(j-1)*m, &m, c+j-1, s+j-1);
      drot_(&ONEi, g+(j-1), &ONEi, g+j, &ONEi, c+j-1, s+j-1);
      mRes[mIter] = fabs( g[j] );

      // Check convergece
      if (mIter == MaxIt-1 || (Type == USUAL && mRes[mIter] < Tol)
          || (Type == REFINE && mRes[mIter] > mRes[mIter-1]/2.0)) {
        QuitIter = true;
        mIter++;
        break;
      }
      
      mIter++;

    }

    // Get the explicit y (stored in g)
    if (QuitIter == false) {
      j = m-1;
    }
    dtrsv_(&UPLO_U, &TRANS_N, &DIAG_N, &j, H, &m, g, &ONEi);
    DVector y(j);
    y.SetBlock(0, j, g);

    // x = x + M(V(:,0:j-1) * y)
    DVector Vy, dx;
    DMatrix VV;
    V.GetBlock(0, N, 0, j, VV);
    VV.MatVec(y, Vy, NORMAL);
    M.MatVec(Vy, dx, NORMAL, mStateM);
    x.Add(dx);

    // Restart
    if (QuitIter) {
      break;
    }
    else {
      // r = b - Ax
      A.MatVec(x, Ax, NORMAL, mStateA);
      b.Subtract(Ax, r);
      // beta = norm(r)
      beta = r.Norm2();
      // reset g
      memset(g, 0, m*sizeof(double));
      g[0] = beta;
    }

  }

  // Clean up
 cleanup:
  Delete_1D_Array<double>(&H);
  Delete_1D_Array<double>(&g);
  Delete_1D_Array<double>(&c);
  Delete_1D_Array<double>(&s);

}


#endif
