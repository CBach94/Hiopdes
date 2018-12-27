#include <armadillo>
using namespace arma;

class Projection{
public:
  Projection(int dim, mat diffmat, vec weights)
  :dim(dim),d(diffmat),w(weights)
  { }
  virtual mat getLaplacian(int dim) = 0;
  virtual mat getDerivativeMatrix(int direction) = 0;
protected:
  int dim;
  mat d;
  vec w;
};

class CollocationProjection : public Projection{
public:
  CollocationProjection(int dim, mat diffmat, vec weights):Projection(dim,diffmat,weights){
  }
  // Returns matrix rep. of [dim] dimensional Laplacian
  mat getLaplacian(){
    return getLaplacian(Projection::dim);
  }
  mat getLaplacian(int dim);
  // Returns matrix rep. of d/dq*\hat{q}, q is a direction
  mat getDerivativeMatrix(int direction);
};

class GalerkinProjection : public Projection{
public:
  GalerkinProjection(int dim, mat diffmat, vec weights):Projection(dim,diffmat,weights){
  }
  // Returns matrix rep. of [dim] dimensional Laplacian
  mat getLaplacian(){
    return getLaplacian(Projection::dim);
  }
  mat getLaplacian(int dim);
  // Returns matrix rep. of d/dq*\hat{q}, q is a direction
  mat getDerivativeMatrix(int direction);
};
