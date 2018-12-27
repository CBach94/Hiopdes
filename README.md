# Hiopdes
HIgh Order Partial Differential Equations Solver framework.
This framework can solve PDEs using galerkin spectral methods.
Currently, following equations can be solves:
- Laplace equation (1D/2D/3D)
- Poisson equation (1D/2D/3D)
- Advection-Diffusion equation (1D/2D/3D)

# Prerequisites
To use the framework you will need
- Armadillo
- Cmake
- GCC

To use the Post-Processing scripts you can use Matlab or Octave.

# Usage
./Solver n N creates a set of files. One file contains the geometrical positions of where the fields were evaluated. Another file contains the fieldvalues at these points.

In octave run for example threeElement n N and the field is plottet on the domain.
![alt text](https://raw.githubusercontent.com/CBach94/Hiopdes/blob/master/Data/AdvDiff.pdf)
