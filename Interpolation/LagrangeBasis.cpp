#include "LagrangeBasis.h"

// Computes L(x) (in 1D)
mat LagrangeBasis::evaluateBasisAt(vec x){
  int n_evalPoints = x.n_rows;

  mat L = zeros(n_evalPoints,n_degree);
  for (size_t j = 0; j < n_degree; j++) {
    vec L_j = ones(n_evalPoints,1);
    for (size_t m = 0; m < n_degree; m++) {
      if (m!=j) {
        L_j%=(x-xq(m))/(xq(j)-xq(m));
      }
    }
    L.col(j) = L_j;
  }

  return L;
}

// Computes dL/dx @ x (in 1D)
mat LagrangeBasis::evaluateBasisDerivativeAt(vec x){
  int n = x.n_rows;
  mat dLdx = zeros(n,n_degree);

  mat d = repmat(xq,1,n_degree) - repmat(xq.t(),n_degree,1);
  // d(i,j) = xq(i) - xq(j);
  for (size_t j = 0; j < n_degree; j++) {
    vec k = zeros(n,1);
    for (size_t l = 0; l < n_degree; l++) {
      if (j!=l) {
        k = ones(n,1)/d(j,l);
        for (size_t m = 0; m < n_degree; m++) {
          if (j!=m && m!=l) {
            k%=(x-xq(m))/d(j,m);
          }
        }
        dLdx.col(j) += k;
      }
    }
  }
  return dLdx.t();
}
