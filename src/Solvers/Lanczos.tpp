#ifndef _LANCZOS_TPP_
#define _LANCZOS_TPP_


//--------------------------------------------------------------------------
template<class MatrixA>
void Lanczos::
Run(MatrixA &A,         // Martix A
    MatState mStateA,   // Needed by A.MatVec()
    const DVector &v,   // Initial vector
    bool PartialReorth, // Whether to perform partial reorth
    INTEGER MaxIt,      // Maximum # of iterations
    bool EarlyStop,     // Whether do early stopping
    INTEGER k,          // Early stopping when k Ritz values converged
    double RTol         // Convergence test for Ritz values
    ) {

  // Norm of v, and residual array
  NormV = v.Norm2();
  Delete_1D_Array<double>(&mRes);
  New_1D_Array<double, INTEGER>(&mRes, MaxIt);

  // Matrices V and T
  INTEGER MaxIt1 = MaxIt + 1;
  INTEGER N = v.GetN();
  V.Init(N, MaxIt1);
  T.Init(MaxIt1);

  // Normalize v and assign to V(:,1)
  DVector w;
  v.Divide(NormV, w);
  V.SetColumn(0, w);

  // Vectors alpha and beta. For computational efficiency and ease of
  // operation, we use the traditional *double array rather than the
  // DVector class.
  double *alpha = NULL, *beta = NULL;
  New_1D_Array<double, INTEGER>(&alpha, MaxIt);
  New_1D_Array<double, INTEGER>(&beta, MaxIt1);
  
  // [If do partial reorth] Matrices G and Omega. For computational
  // efficiency and ease of operation, we use the traditional *double
  // array rather than the DMatrix class.
  //
  // Note: The matrices are transposed from the matlab code to enhance
  // locality.
  double *G = NULL, *Omega = NULL;
  if (PartialReorth) {
    New_1D_Array<double, INTEGER>(&G, 4 * MaxIt);
    New_1D_Array<double, INTEGER>(&Omega, MaxIt1 * MaxIt1);
  }

  // Needed by reorth
  double Cons = sqrt((double)N) * EPS_2;
  double PartialReorthTol = SqrtEPS;
  bool NeedPartialReorth = false;

  // Iteration count
  mIter = 0;
  bool QuitIter = false;

  // Loop
  DVector vi, vi1;
  DMatrix V2, T2;
  for (INTEGER i = 0; i < MaxIt; i++) {

    // w = A(V(:,i)) - beta(i)*V(:,i-1)
    V.GetColumn(i, vi);
    A.MatVec(vi, w, NORMAL, mStateA);
    if (i > 0) {
      DVector bvi1;
      vi1.Multiply(beta[i], bvi1);
      w.Subtract(bvi1);
    }

    // alpha(i) = w' * V(:,i)
    alpha[i] = w.InProd(vi);

    // w = w - alpha(i)*V(:,i)
    DVector avi;
    vi.Multiply(alpha[i], avi);
    w.Subtract(avi);

    // Reorthogonalization
    if (PartialReorth) {

      if (i == 0) {
        G[1] = 0.0;
        G[2] = 0.0;
        G[3] = alpha[0] * alpha[0];
        G[0] = G[3];
        Omega[0] = 1.0;
        Omega[0+MaxIt1] = Cons;
        Omega[1+MaxIt1] = 1.0;
      }
      else {
        beta[i+1] = w.Norm2();
        G[1+i*4] = G[1+i*4] + fabs( beta[i-1] * beta[i] );
        G[2+i*4] = G[2+(i-1)*4] + beta[i] * beta[i] +
          fabs( (alpha[i-1] + alpha[i]) * beta[i] );
        G[3+i*4] = beta[i] * beta[i] + alpha[i] * alpha[i] +
          fabs( (alpha[i-1] + alpha[i]) * beta[i] ) +
          fabs( beta[i-1] * beta[i] );
        G[0+i*4] = Max<double>(Max<double>(Max<double>(G[0+(i-1)*4], G[1+i*4]),
                                           G[2+i*4]), G[3+i*4]);
        double theta = Cons * sqrt(G[0+i*4]);
        INTEGER j = 0;
        Omega[j+(i+1)*MaxIt1] = ( beta[j+1] * Omega[(j+1)+i*MaxIt1] +
                                  (alpha[j] - alpha[i]) * Omega[j+i*MaxIt1] -
                                  beta[i] * Omega[j+(i-1)*MaxIt1] + theta ) /
                                beta[i+1];
        for (j = 1; j <= i-1; j++) {
          Omega[j+(i+1)*MaxIt1] = ( beta[j+1] * Omega[(j+1)+i*MaxIt1] +
                                    (alpha[j] - alpha[i]) * Omega[j+i*MaxIt1] +
                                    beta[j] * Omega[(j-1)+i*MaxIt1] -
                                    beta[i] * Omega[j+(i-1)*MaxIt1] + theta ) /
                                  beta[i+1];
        }
        Omega[i+(i+1)*MaxIt1] = Cons;
        Omega[(i+1)+(i+1)*MaxIt1] = 1.0;
      }

      INTEGER j = 0;
      double MaxAbsOmega = fabs(Omega[j+(i+1)*MaxIt1]);
      for (j = 1; j <= i; j++) {
        double tm = fabs(Omega[j+(i+1)*MaxIt1]);
        MaxAbsOmega = MaxAbsOmega > tm ? MaxAbsOmega : tm;
      }
      if ( (i>0) && (NeedPartialReorth || MaxAbsOmega >= PartialReorthTol) ) {
        // Reorth w = w - V(:,0:i-1) * (V(:,0:i-1)' * w)
        V.GetBlock(0, N, 0, i, V2);
        DVector w1, w2;
        V2.MatVec(w, w1, TRANSPOSE);
        V2.MatVec(w1, w2, NORMAL);
        w.Subtract(w2);
        // Book-keeping
        for (j = 0; j <= i; j++) {
          Omega[j+(i+1)*MaxIt1] = Cons;
        }
        Omega[(i+1)+(i+1)*MaxIt1] = 1.0;
        if (NeedPartialReorth == true) {
          NeedPartialReorth = false;
        }
        else {
          NeedPartialReorth = true;
        }
      }

    }
    else {

      // Reorth w = w - V(:,0:i-1) * (V(:,0:i-1)' * w)
      if (i > 0) {
        V.GetBlock(0, N, 0, i, V2);
        DVector w1, w2;
        V2.MatVec(w, w1, TRANSPOSE);
        V2.MatVec(w1, w2, NORMAL);
        w.Subtract(w2);
      }

    }

    // beta(i+1) = norm(w)
    beta[i+1] = w.Norm2();

    // V(:,i+1) = w/beta(i+1)
    w.Divide(beta[i+1]);
    V.SetColumn(i+1, w);

    // T(i,i) = alpha(i)
    // T(i+1,i) = beta(i+1)
    // T(i,i+1) = beta(i+1)
    T.SetEntry(i, i, alpha[i]);
    T.SetEntry(i+1, i, beta[i+1]);
    T.SetEntry(i, i+1, beta[i+1]);

    // Book-keeping
    mIter++;

    // Check convergence of Ritz values (do not check in every step)
    if ((EarlyStop && (mIter%k == 0)) || (mIter == MaxIt)) {

      // T(0:i,0:i) = W * S * W'
      T.GetBlock(0, mIter, 0, mIter, T2);
      DMatrix W;
      T2.SymEig(S, W);

      // Res(j) = abs( beta(i+1) * W(i,j) ). Also check convergence
      QuitIter = true;
      for (INTEGER j = 0; j < mIter; j++) {
        mRes[j] = fabs( beta[i+1] * W.GetEntry(i,j) );
        if (j >= mIter - k && mRes[j] > RTol) {
          // Eigenvalues are ordered increasingly
          QuitIter = false;
        }
      }

      // If converged, compute Ritz vectors U = V(:,0:i) * W and exit
      if (QuitIter || mIter == MaxIt) {
        V.GetBlock(0, N, 0, mIter, V2);
        V2.MatMat(W, U, NORMAL, NORMAL);
        break;
      }

    }

    // Book-keeping
    vi1 = vi;

  }

  // On output
  //
  // Trim matrices
  V.GetBlock(0, N, 0, mIter, V2);
  V = V2;
  T.GetBlock(0, mIter, 0, mIter, T2);
  T = T2;

  // Clean up
  Delete_1D_Array<double>(&alpha);
  Delete_1D_Array<double>(&beta);
  Delete_1D_Array<double>(&G);
  Delete_1D_Array<double>(&Omega);

}


#endif
