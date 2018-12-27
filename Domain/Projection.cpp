#include "Projection.h"

mat CollocationProjection::getLaplacian(int dim){
  mat d2 = d*d;
  if (dim==1) {
    return d2;
  } else{
    mat L = d2;
    mat I = eye(L.n_rows,L.n_rows);
    for (size_t i = 0; i < dim-1; i++) {
      L = kron(L,I);
    }
    return L + kron(I,getLaplacian(dim-1));
  }
}

mat CollocationProjection::getDerivativeMatrix(int direction){
  if (dim==1){
    return d;
  } else {
    mat I = eye(d.n_rows,d.n_rows);
    mat diffmat = I;
    for (size_t i = 0; i < dim-1; i++) {
      if (i==0) {
        if (i==direction) {
          diffmat = kron(d,diffmat);
        } else {
          diffmat = kron(diffmat,d);
        }
      } else{
        diffmat = kron(diffmat,I);
      }
    }
  }
}

mat GalerkinProjection::getLaplacian(int dim){
  mat W = diagmat(w);
  mat d2 = d.t()*W*d;
  if (dim==1) {// if dimension == 1, derivative is done wrt. the only coord.
    return d2;
  } else{// Otherwise recursively construct the matrix rep. of the operator
    mat L = d2;
    mat W = diagmat(w);
    for (size_t i = 0; i < dim-1; i++) {
      L = kron(L,W);
    }
    return -1.0*(L + kron(W,getLaplacian(dim-1)));
  }
}
mat GalerkinProjection::getDerivativeMatrix(int direction){
  mat W = diagmat(w);

  if (dim==1){
    return -d.t()*W;
  } else {
    mat diffmat = W;
    mat D = d.t()*W;
    for (size_t i = 0; i < dim-1; i++) {
      if (i==0) {
        if (i==direction) {
          diffmat = kron(D,diffmat);
        } else {
          diffmat = kron(diffmat,D);
        }
      } else{
        diffmat = kron(diffmat,W);
      }
    }
    return -diffmat;
  }
}
