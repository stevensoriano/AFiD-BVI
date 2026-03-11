# AFiD-BVI

A highly parallel application for simulating Blade-Vortex Interactions (BVI) using direct numerical simulation (DNS).

AFiD-BVI is a code for parallel numerical simulations of Blade-Vortex Interactions (BVI) built from [AFiD](https://github.com/PhysicsofFluids/AFiD). The code solves the incompressible Navier-Stokes equations for a Newtonian fluid coupled with an immersed body through the Immersed Boundary Method (IBM). It was developed to study the interaction of a thin cylinder (wire) cutting through a columnar vortex tube at normal incidence, as described in [Soriano & Ostilla-Mónico (2024)](#citing-this-work).

## Physical Problem

The code simulates normal blade-vortex interaction (BVI): a cylinder of diameter $D$ translates at velocity $V$ through a Lamb-Oseen vortex of circulation $\Gamma$ and core radius $\sigma_0$. The vortex is oriented perpendicular to the cylinder axis. The interaction is governed by three dimensionless parameters:

- **Impact parameter**: $I_P = 2\pi\sigma_0 V / \Gamma$ — ratio of the body's free-stream velocity to the maximum vortex swirl velocity
- **Circulation Reynolds number**: $Re_\Gamma = \Gamma / \nu$ — set via `REN` in `bou.in`
- **Thickness parameter**: $T = D / \sigma_0$ — ratio of the cylinder diameter to the vortex core radius

Two interaction regimes are distinguished:

| Regime | Condition | Behavior |
|--------|-----------|----------|
| **Strong vortex** | $I_P < 0.08$ | Boundary layer separates and wraps around the primary vortex before impact, heavily disrupting the vortex |
| **Weak vortex** | $I_P > 0.2$ | Minimal boundary layer interaction; the vortex deforms and propagates axial waves but retains coherence |

## Numerical Methods

The code uses energy-conserving second-order centered finite differences for spatial discretization. Time marching is done with third-order Runge-Kutta for the nonlinear terms and second-order Crank-Nicolson for the viscous terms. The Immersed Boundary Method (IBM) with a moving least squares (MLS) formulation is used to incorporate the body into the flow field.

The computational domain is a triply periodic box. The body geometry is specified in GTS (GNU Triangulated Surface) format. The vortex is initialized using a Lamb-Oseen velocity profile with circulation $\Gamma = 1$ (hardcoded in `CreateInitialConditions.F90`) and core radius $\sigma_0$ set via `VORTR` in `bou.in`. The cylinder velocity $V$ is set via `VORTSPEED` in `bou.in` (see the [manual](MANUAL.md) for details on parameter naming).

## Features

- Provision for handling multiple deformable bodies represented using an Immersed Boundary approach
- An interaction-potential approach allows for the deformation of the immersed bodies
- Moving Least Squares (MLS) for information exchange between the fluid and immersed body
- Pencil-type domain decomposition
- Hybrid MPI/OpenMP parallelization
- [2DECOMP&FFT](https://github.com/2decomp-fft/2decomp-fft) library for global data transpositions

## Usage

Consult the [manual](MANUAL.md) for details on the code input parameters and their relationship to the dimensionless groups in the paper. An example usage of the code can be found in [`example/AFiD_BVI_Example.pdf`](https://github.com/stevensoriano/AFiD-BVI/blob/main/example/AFiD_BVI_Example.pdf).

The default configuration in `bou.in` corresponds to the $I_P = 0.25$, $Re_\Gamma = 1000$ case (weak vortex regime) from Table I of Soriano & Ostilla-Mónico (2024).

## Installation

### Prerequisites

* MPI
* BLAS
* LAPACK
* FFTW3
* HDF5 with parallel I/O

The Makefile on the main page was utilized for code compilation on the Carya system, part of the Research Computing Data Core (RCDC) at the University of Houston. Prior to compilation, the following specific modules were loaded.

```bash
module load intel
module load FFTW
module load HDF5
```

## Citing This Work

If you use this code, please cite the following:

**Primary paper:**

Soriano, S. & Ostilla-Mónico, R. "Direct numerical simulations of a cylinder cutting a vortex." *Phys. Rev. Fluids*, **9**, 054701 (2024). https://doi.org/10.1103/PhysRevFluids.9.054701

```bibtex
@article{PhysRevFluids.9.054701,
  title   = {Direct numerical simulations of a cylinder cutting a vortex},
  author  = {Soriano, Steven and Ostilla-M\'onico, Rodolfo},
  journal = {Phys. Rev. Fluids},
  volume  = {9},
  issue   = {5},
  pages   = {054701},
  year    = {2024},
  month   = {May},
  publisher = {American Physical Society},
  doi     = {10.1103/PhysRevFluids.9.054701},
  url     = {https://link.aps.org/doi/10.1103/PhysRevFluids.9.054701}
}
```

**Base numerical code (AFiD):**

Van Der Poel, E. P., Ostilla-Mónico, R., Donners, J. & Verzicco, R. "A pencil distributed finite difference code for strongly turbulent wall-bounded flows." *Computers & Fluids*, **116**, 10-16 (2015). https://doi.org/10.1016/j.compfluid.2015.04.007

**Immersed boundary method (MLS formulation):**

Spandan, V., Meschini, V., Ostilla-Mónico, R., Lohse, D., Querzoli, G., de Tullio, M. D. & Verzicco, R. "A parallel interaction potential approach coupled with the immersed boundary method for fully resolved simulations of deformable interfaces and membranes." *Journal of Computational Physics*, **348**, 567-590 (2017). https://doi.org/10.1016/j.jcp.2017.07.036
