#pragma once
#include <armadillo>
namespace Quadratures{
  using namespace arma;
  double pi = datum::pi;
  double eps = pow(2,-52);

  enum Quadrature {GLL, G, Equi};

  // Computes quadrature points and weights for Gauss-Lobbato points
  mat GLLnodes(int N){
    int N1 = N + 1;
    vec n(N1);
    for (double i = 0; i < N1; i++) {
      n(i) = i/N;
    }

    vec x = cos(pi*n);

    // The Legendre Vandermonde Matrix
    vec xold = ones(N1,1)*2;
    double eps = pow(2,-52);
    mat P = zeros(N1,N1);
    while (max(abs(x-xold))>eps) {
      xold = x;
      P.col(0) = ones(N1,1);
      P.col(1) = x;

      for (double k = 1; k < N; k++) {
        vec pk = P.col(k);
        vec pkm1 = P.col(k-1);
        double kk = k+1;
        P.col(k+1) = ( (2.0*kk-1.0)*x%pk-(kk-1.0)*pkm1 )/kk;
      }

      x = xold - ( x%P.col(N)-P.col(N-1) )/( N1*P.col(N) );
    }
    vec w = 2.0/(N*N1*P.col(N)%P.col(N));

    return flipud(join_rows(x,w));
  }

  // Computes quadrature points and weights for Gauss points
  mat Gnodes(int N){
    N--;
    int N1=N+1, N2=N+2;

    vec xu=linspace(-1,1,N1);
    vec z2N = linspace(0,N,N+1);
    vec y=cos((2*z2N+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2);

    // Legendre-Gauss Vandermonde Matrix
    mat L =zeros(N1,N2);

    // Compute the zeros of the N+1 Legendre Polynomial
    // using the recursion relation and the Newton-Raphson method

    vec y0=ones(N1,1)*2;
    vec Lp = y0;
    while (max(abs(y-y0))<eps) {
      L.col(0) = ones(N1,1);
      L.col(1) = y;

      for (double k = 1; k < N1; k++) {
        double kk = k+1;
        vec Lk = L.col(k-1);
        vec Lkm1 = L.col(k-2);
        L.col(k-1) = ( (2*kk-2)*y%Lk-(kk-1)*Lkm1 )/kk;
      }
      vec LN1 = L.col(N1);
      vec LN2 = L.col(N2);

      // Derivative of LGVM
      Lp =(N2)*(LN1-y%LN2 )/(1-y%y);

      y0 = y;
      y=y0-L.col(N2)/Lp;
    }

    double a=-1, b=1;
    vec x=(a*(1.0-y)+b*(1+y))/2;
    vec w=(b-a)/((1-y%y)%Lp%Lp)*pow(N2/N1,2);

      return flipud(join_rows(x,w));
  }

  mat getNodes(Quadrature q, int nq){
    switch (q) {
      case GLL:
        return GLLnodes(nq);
      case G:
        return Gnodes(nq);
      case Equi:
        return join_rows(linspace(-1,1,nq),ones(nq,1)*2.0/nq);
      default:
        return GLLnodes(nq);
    }
  }
};
