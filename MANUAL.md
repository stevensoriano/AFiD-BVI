# AFiD-BVI Manual

This manual documents the input parameters, code structure, and output files for AFiD-BVI. For the physical background and simulation results, see [Soriano & Ostilla-Mónico, Phys. Rev. Fluids 9, 054701 (2024)](https://doi.org/10.1103/PhysRevFluids.9.054701).

## Relationship Between Input Parameters and Dimensionless Groups

The vortex circulation is hardcoded as &Gamma; = 1.0 in `CreateInitialConditions.F90`. The key dimensionless groups from the paper are recovered from the input parameters as follows:

| Dimensionless group | Formula | Input parameters |
|---------------------|---------|------------------|
| Impact parameter I_P | 2&pi;&sigma;_0 V / &Gamma; | `VORTR` (&sigma;_0), `VORTSPEED` (V), &Gamma; = 1 |
| Circulation Reynolds number Re_&Gamma; | &Gamma; / &nu; | `REN` = Re_&Gamma; (since &Gamma; = 1, &nu; = 1/REN) |
| Thickness parameter T | D / &sigma;_0 | D from GTS geometry scaled by `SCLF`, `VORTR` (&sigma;_0) |
| Cylinder Reynolds number Re_c | V D / &nu; | Re_c = Re_&Gamma; I_P T / (2&pi;) |

**Example (default `bou.in`):** With `VORTR` = 0.095, `VORTSPEED` = 0.419, `REN` = 1000, and &Gamma; = 1:
- I_P = 2&pi; &times; 0.095 &times; 0.419 / 1.0 = 0.25
- Re_&Gamma; = 1000
- This corresponds to the weak vortex regime case (I_P = 0.25, Re_&Gamma; = 1000) from Table I of the paper.

To simulate other cases from the paper, adjust `VORTSPEED` (V) and `REN` (Re_&Gamma;). For example, for I_P = 0.05 at Re_&Gamma; = 1000: V = I_P &times; &Gamma; / (2&pi;&sigma;_0) = 0.05 / (2&pi; &times; 0.095) = 0.0838.

## Code Structure

The main loop with the integration scheme can be found in `TimeMarcher.F90`. It calls the following routines:

 * `ExplicitTermsVX`
 * `ExplicitTermsVY`
 * `ExplicitTermsVZ`
 * `ExplicitTermsTemp`
 * `ImplicitAndUpdateVX`
 * `ImplicitAndUpdateVY`
 * `ImplicitAndUpdateVZ`
 * `ImplicitAndUpdateTemp`
 * `update_halo`
 * `CalculateLocalDivergence`
 * `SolvePressureCorrection`
 * `CorrectVelocity`
 * `CorrectPressure`
 * `CalcMLSForce`
 * `MLSForceVelocity`

## List of Input and Output Files

### Input files (.in)
 * `bou.in` — Main simulation parameters (see next section)
 * `mlspart.in` — Input variables for the immersed bodies (see next section)
 * `spos.in` — Initial position (x, y, z) of the body/bodies in the computational domain
 * `ppart.in` — Point particle input variables
 * `stst2.in` — Statistics probe locations and slab dump positions

### Output data directories
 * `flowmov/` — HDF5 files of flow field data (vorticity, velocity, and pressure)
 * `vtkfiles/` — VTK files of body force data

### Output data files (.h5)
 * `continua_vx.h5`
 * `continua_vy.h5`
 * `continua_vz.h5`
 * `continua_temp.h5`
 * `continua_master.h5`
 * `cordin_info.h5`

### Output log files (.out)
 * `Total_time.out`
 * `dissipation.out` — Numerical dissipation (used to check resolution adequacy)
 * `fforce.out` — Force on the immersed body
 * `fforce_nor.out` — Normalized force on the immersed body
 * `rms_vel.out` — RMS velocity statistics
 * `viscf.out` — Viscous force data

## Input Variables in `bou.in`

**Note on naming:** The variable labels in `bou.in` (`VORTR`, `VORTSPEED`) do not match the internal Fortran variable names. In the source code (`ReadInputFile.F90`), these are read as `vortexrad` and `wirev` respectively.

```
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
bou.in
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
N1M             : Number of computational cells in the x-direction
N2M             : Number of computational cells in the y-direction
N3M             : Number of computational cells in the z-direction
NSST            : if (NSST = 3) third-order Runge-Kutta integrator
                  if (NSST = 1) Adams-Bashforth integrator
NREAD           : if (NREAD = 1) Start simulation by reading continuation files
                  if (NREAD = 0) Start from initial conditions created in the code

NTST            : Total number of time steps
WALLTIMEMAX     : Maximum wall time of the simulation (seconds)
IRESET          : if (IRESET = 1) Simulation time resets to 0
                  if (IRESET = 0) Simulation time unchanged

XLEN, YLEN, ZLEN: Lengths of the periodic domain in each direction.
                   The domain corresponds to approximately
                   (XLEN/VORTR) x (YLEN/VORTR) x (ZLEN/VORTR) vortex radii.
                   Default: 3.0 x 3.0 x 9.0, giving ~31.6 sigma_0 x 31.6 sigma_0 x 94.7 sigma_0

REN             : Circulation-based Reynolds number, Re_Gamma = Gamma / nu.
                  Since Gamma = 1.0 (hardcoded), nu = 1 / REN.
DT              : Initial time step
RESID           : Maximum allowable residual for mass conservation
CFLMAX          : Maximum CFL number for simulation (recommended: 1.2)

TSTA            : Time to start collecting statistics
STAREAD         : if (STAREAD = 1) Read statistics from continuation files
                  if (STAREAD = 0) Don't read statistics

IDTV            : if (IDTV = 1) Variable time stepping
                  if (IDTV = 0) Fixed time stepping
DTMAX           : Maximum time step for variable time stepping

MOVIE           : if (MOVIE = 1) Flow field is saved for visualization (HDF5 files in flowmov/)
                  if (MOVIE = 0) Flow field is not saved
DUMPMOVIET      : Time interval between flow field snapshots

VORTR           : Vortex core radius, sigma_0.
                  This is the Lamb-Oseen core radius used in CreateInitialConditions.F90.
                  Internal variable name: vortexrad
VORTSPEED       : Cylinder (wire) translational velocity, V.
                  Despite the name, this is NOT the vortex speed or circulation.
                  The impact parameter is computed as I_P = 2*pi*VORTR*VORTSPEED / Gamma,
                  where Gamma = 1.0 is hardcoded.
                  Internal variable name: wirev

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
```

## Input Variables in `mlspart.in`

```
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
mlspart.in
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
IMLSFORC        : if (IMLSFORC = 1) MLS forcing onto the flow is turned on
                  if (IMLSFORC = 0) MLS forcing onto the flow is turned off
IMLSSTRUC       : if (IMLSSTRUC = 1) MLS forcing onto the object is turned on
                  if (IMLSSTRUC = 0) MLS forcing onto the object is turned off
PREAD           : if (PREAD = 1) Read body positions from restart file
                  if (PREAD = 0) Position bodies at spos.in coordinates

SCLF            : Scale factor for the GTS geometry.
                  The cylinder diameter D in simulation units is determined by
                  the GTS file geometry multiplied by this scale factor.
                  The thickness parameter is T = D / VORTR.
RHOP            : Density ratio of the body to the fluid

GTSFX           : GTS geometry file to be read in (e.g., cylinder.gts)
DATFX           : Name of restart file for the object positions

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
```
