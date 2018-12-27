#include "Quadratures.h"
#include "Projection.h"

using namespace Quadratures;

template <int dimension, class Interpolation>
class Element{
public:
  Element(Quadrature quadType, int n_degree)
  :interpolator(Interpolation(getNodes(quadType,n_degree).col(0)))
  {
    mat xw = getNodes(quadType,n_degree);
    quadPoints  = xw.col(0);
    quadWeights = xw.col(1);

  }

  // Imposes a Dirichlet condition on an equation system Ax=b
  void imposeDirichlet(mat & A,vec & b, uvec indicies, double val){
    A.rows(indicies).fill(0);
    b(indicies).fill(val);
    for (size_t i = 0; i < indicies.n_rows; i++) {
      A(indicies(i),indicies(i)) = 1;
    }
  }
  void imposeDirichlet(mat & A,vec & b, uvec indicies, vec bVals){
    A.rows(indicies).fill(0);
    b(indicies)=(bVals);
    for (size_t i = 0; i < indicies.n_rows; i++) {
      A(indicies(i),indicies(i)) = 1;
    }
  }

protected:
  Interpolation interpolator;
  vec quadPoints;
  vec quadWeights;
};

enum Region{North,South,East,West,Front,Back};

template<class Interpolation, class Projection>
class Std2DElement : public Element<2,Interpolation>{
public:
  using ThisElmt = Element<2,Interpolation>;
  Std2DElement(Quadrature quadType, int n_degree)
  :Element<2,Interpolation>(quadType,n_degree),
  quadType(quadType),
  n_degree(n_degree)
  {
    int n = n_degree;
    // Boundary indicies;
    uvec iNorth = regspace<uvec>(0,n);
    uvec iSouth = iNorth + n*(n+1);
    uvec iWest  = regspace<uvec>(0,n+1,n*(n+1));
    uvec iEast  = iWest + n;

    boundaryIndicies[0] = iEast;
    boundaryIndicies[1] = iWest;
    boundaryIndicies[2] = iNorth;
    boundaryIndicies[3] = iSouth;
  }
  mat getLaplacian(){
    mat xw = getNodes(quadType,n_degree);
    auto projector = Projection(
      2, // dimensino
      ThisElmt::interpolator.getDiffMat(),
      xw.col(1)
    );
    return projector.getLaplacian();
  }

  mat getDerivativeMatrix(int direction){
    mat xw = getNodes(quadType,n_degree);
    auto projector = Projection(
      2, // dimensino
      ThisElmt::interpolator.getDiffMat(),
      xw.col(1)
    );
    return projector.getDerivativeMatrix(direction);
  }

  // Create a global matrix system from two elements
  void connectElements(
    mat & LHSi,      // matrix element i
    vec & RHSi,      // source element i
    uvec indicies_i, // indicies in element i that Connects with element j
    mat & LHSj,      // matrix element j
    vec & RHSj,      // source element j
    uvec indicies_j  // indicies in element j that Connects with element i
  ){
    // Create global matrix
      int ni = LHSi.n_rows;
      int nj = LHSj.n_rows;
      mat A = zeros(ni+nj,ni+nj);

    // off-diagonal block matrices
      mat sharedi = zeros(ni,nj);
      mat sharedj = zeros(nj,ni);

      sharedi.rows(indicies_i) = LHSj.rows(indicies_j);
      sharedj.rows(indicies_j) = LHSi.rows(indicies_i);

      A.rows(0,ni-1)     = join_rows(LHSi,sharedi);
      A.rows(ni,ni+nj-1) = join_rows(sharedj,LHSj);

      int nI = indicies_i.n_rows;
      for (size_t i = 0; i < nI; i++) {
        int Ii = indicies_i(i);
        int Ij = ni + indicies_j(i);

        A.row(Ii).fill(0);

        A(Ii,Ii) = 1;
        A(Ii,Ij) = -1;
      }
    // Global source term
      vec b = join_cols(RHSi,RHSj);

    LHSi = A;
    RHSi = b;
  }
  void setDirichletAt(uvec indicies, mat & A, vec & b, vec bVals){
    ThisElmt::imposeDirichlet(indicies, A, b, bVals);
  }
  void setDirichletAt(Region region, mat & A, vec & b, double val){
    ThisElmt::imposeDirichlet(A,b,getIndiciesAt(region),val);
  }
  uvec getIndiciesAt(Region region){
    switch (region) {
      case North:
        return boundaryIndicies[0];
      case South:
        return boundaryIndicies[1];
      case West:
        return boundaryIndicies[2];
      case East:
        return boundaryIndicies[3];
    }
  }
  vec evaluateFieldAt(vec field, vec points){
    mat L = ThisElmt::interpolator.evaluateAt(points);
    mat H = kron(L,L);

    return H*field;
  }
private:
  uvec boundaryIndicies[4];
  Quadrature quadType;
  int n_degree;
};
