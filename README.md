# AFiD-BVI
A highly parallel application for Blade-Vortex Interactions (BVI). 

AFiD-BVI is a code for parallel numerical simulations of Blade-Vortex Interactions (BVI) built from https://github.com/PhysicsofFluids/AFiD. The code will solve the fluid flow of an incompressible Newtonian fluid coupled with fluid-structure interactions.

The code uses energy-conserving centered finite differences to discretize the domain spatially. Time marching is done third-order Runge-Kutta for the non-linear terms and second-order Crank-Nicolson for the viscous terms. Using a moving least squares formulation, the Immersed Boundary Method (IBM) is used to incorporate the body into the flow field.

The computational domain is a triply periodic box that allows a body of any geometry to be incorporated into the flow field. The body is represented in a GTS (GNU Triangulated Surface) file format. The body can be programmed to move within the flow field. The Vortex is initialized using a Lamb–Oseen vortex distribution and can be oriented to solve normal, parallel, and oblique BVI. 

**Reference**

Van Der Poel, Erwin P., et al. "A pencil distributed finite difference code for strongly turbulent wall-bounded flows." *Computers & Fluids* 116 (2015): 10-16. https://www.sciencedirect.com/science/article/pii/S0045793015001164

Spandan, Vamsi, et al. "A parallel interaction potential approach coupled with the immersed boundary method for fully resolved simulations of deformable interfaces and membranes." *Journal of Computational Physics* 348 (2017): 567-590. https://www.sciencedirect.com/science/article/pii/S0021999117305442

## Usage

Consult the manual for details on the code input parameters. An example usage of the code can be found here.

## Installation 

### Prerequisites

* MPI
* BLAS
* LAPACK
* FFTW3
* HDF5 with parallel I/O

The provided Makefile was used to compile the code on Carya at the Research Computing Data Core (RCDC) at the University of Houston. The following modules were loaded before compiling.

```bash
module load intel
module load FFTW
module load HDF5
```




 



