#include <armadillo>

using namespace arma;
class LagrangeBasis{
public:
  LagrangeBasis(vec x):xq(x)
  {
    n_degree = xq.n_rows;

    // only true on quadrature points due to cardinality
    basisFunctions           = eye(n_degree,n_degree);
    // compute derivatives
    basisFunctionDerivatives = evaluateBasisDerivativeAt(xq);
  }
  int getDegree(){return n_degree;}
  vec getQuadraturePoints(){return xq;}
  mat getBasisFunctionDerivatives(){ return basisFunctionDerivatives;}
  mat evaluateBasisAt(vec x);
  mat evaluateBasisDerivativeAt(vec x);
private:
  vec xq;
  mat basisFunctions;
  mat basisFunctionDerivatives;
  int n_degree;
};
