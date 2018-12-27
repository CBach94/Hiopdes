#include <armadillo>

using namespace arma;

template<class InterpolationBasis>
class Interpolation{
public:
  Interpolation(vec x):interpolator(InterpolationBasis(x))
  { }

  mat evaluateAt(vec x){ return interpolator.evaluateBasisAt(x); }
  mat evaluateDerivativeAt(vec x){
    return interpolator.evaluateBasisDerivativeAt(x);
  }
  mat getDiffMat(){ return interpolator.getBasisFunctionDerivatives().t();}
  vec getQuadraturePoints(){return interpolator.getQuadraturePoints();}
private:
  InterpolationBasis interpolator;
};
