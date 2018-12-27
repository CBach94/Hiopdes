#include "Interpolation/Interpolation.h"
#include "Interpolation/LagrangeBasis.h"
#include "Domain/Element.h"
#include <stdlib.h>

int main(int argc, char const *argv[]) {
  using namespace arma;

  // n: no. quadrature points ~ input argument 1
  // N: no. evaluation points ~ input argument 2
  int n, N;
  if (argc>=3) {
    n = strtol(argv[1], NULL, 10);
    N = strtol(argv[2], NULL, 10);
  } else {
    n = 10;
    N = 10;
  }
  wall_clock timer;
  cout << "n=" << n << "\t N="<<N<< endl;

  using interpolationMethod = Interpolation<LagrangeBasis>;
  using projectionMethod    = GalerkinProjection;

  // Timer for matrix overhead and solving.
  timer.tic();
  // Create a 2D standard element ~ [-1,1]^2
  auto element = Std2DElement<interpolationMethod,projectionMethod>(GLL,n);

  // Create LHS matrix
    // Diffusion term
    double mu = 0.05; // Diffusion constant
    mat Laplacian = element.getLaplacian();

    // Advective term
    vec u = {1.0, 1.0}; // Advective velocity
    mat Advection = u(0)*element.getDerivativeMatrix(0)
                  + u(1)*element.getDerivativeMatrix(1);

  // Alloc matrices and vectors for problems
    // Single element problem
      // PDE operator
      mat LHS = mu*Laplacian - Advection;
      // Source term
      vec RHS = zeros((n+1)*(n+1),1);

    // Three element problem
      // PDE operator
      mat LHS_top    = LHS; // top element
      mat LHS_bottom = LHS; // bottom elm.
      mat LHS_topR   = LHS; // top-right elm.
      // Source term
      vec RHS_top    = RHS; // top elm.
      vec RHS_topR   = RHS; // top-right elm.
      vec RHS_bottom = RHS; // bottom elm.

  // Impose boundary conditions
    element.setDirichletAt(North,LHS_top ,RHS_top,100);
    element.setDirichletAt(West ,LHS_top ,RHS_top,100);

    element.setDirichletAt(South, LHS_bottom, RHS_bottom,0  );
    element.setDirichletAt(West , LHS_bottom, RHS_bottom,100);
    element.setDirichletAt(East , LHS_bottom, RHS_bottom,0  );

    element.setDirichletAt(North ,LHS_topR ,RHS_topR, 100);
    element.setDirichletAt(South ,LHS_topR ,RHS_topR, 0  );
    element.setDirichletAt(East  ,LHS_topR ,RHS_topR, 0  );

    // Connect elements
    element.connectElements(
      LHS_top,    RHS_top,    element.getIndiciesAt(South),
      LHS_bottom, RHS_bottom, element.getIndiciesAt(North)
    );
    element.connectElements(
      LHS_top,  RHS_top,  element.getIndiciesAt(East),
      LHS_topR, RHS_topR, element.getIndiciesAt(West)
    );

    // Matrix system recast
    mat A = LHS_top;
    vec b = RHS_top;

  // Solve PDEs
    vec psi = solve(A,b);
    double time = timer.toc();
    cout << "Elapsed time: " << time << endl;

  // Output/Post-Processing
    mat Psi = reshape(psi,psi.n_rows/3,3);

    // GLL spaced points example
    mat x = getNodes(GLL,n).col(0);
    x.save("xe.txt", raw_ascii);

    vec psiT  = Psi.col(0);
    vec psiB  = Psi.col(1);
    vec psiTR = Psi.col(2);

    // Two element problem
    psiB.save(  "psiB.txt"  ,raw_ascii);
    psiT.save(  "psiT.txt"  ,raw_ascii);
    psiTR.save( "psiTR.txt" ,raw_ascii);

  return 0;
}
