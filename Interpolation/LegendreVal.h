#include<armadillo>
/*
  Evaluates the Legendre polynomial of order N and its derivative at the points in vector x.
  Returns a matrix with the first column being the polynomial and the second column being its derivative
*/
using namespace arma;
mat LegendreVal(vec x,int N){
  int nx = x.n_rows;

  mat L = zeros(nx,N+1);
  L.col(0) = ones(nx,1);

  mat dL = zeros(nx,N+1);
  if(N>0){
    L.col(1) = x;
    dL.col(1) = ones(nx,1);
  }

  for (double k = 1; k < N; k++) {
    double kk= k+1;
    vec Lk = L.col(k);
    vec Lkm1 = L.col(k-1);

    L.col(k+1)  = (2*kk-1)/kk*x%Lk-(kk-1)/kk*Lkm1;

    vec dLkm1 = dL.col(k-1);
    dL.col(k+1) = dLkm1 + (2*kk-1)*Lk;
  }
  return join_rows(L.col(N),dL.col(N));
}
