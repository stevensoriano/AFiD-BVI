Code structure
--------------

The main loop with the integration scheme can be found in TimeMarcher.F90. It calls the following routines:

 * ExplicitTermsVX
 * ExplicitTermsVY
 * ExplicitTermsVZ
 * ExplicitTermsTemp
 * ImplicitAndUpdateVX
 * ImplicitAndUpdateVY
 * ImplicitAndUpdateVZ
 * ImplicitAndUpdateTemp
 * update_halo
 * CalculateLocalDivergence
 * SolvePressureCorrection
 * CorrectVelocity
 * CorrectPressure
 * CalcMLSForce
 * MLSForceVelocity

List of input and output files
------------------------------

Input files (.in)
 *  bou.in - Detailed information on the input variables is given in the next section
 *  mlspart.in - Input variables for the immersed bodies, given in the next section
 *  spos.in - Initial position of the body/bodies in the computational domain
 *  ppart.in - Point particle input variables
 *  stst2.in 

Output data directories
* flowmov - directory containing h5 files of flow field data (vorticity, velocity, and pressure)
* vtkfiles - directory containing vtk files of body force data 
  
Output data files (.h5)
 * continua_vx.h5
 * continua_vy.h5
 * continua_vz.h5
 * continua_temp.h5
 * continua_master.h5
 * cordin_info.h5 

Output log files (.out)
 * Total_time.out
 * dissipation.out
 * fforce.out
 * fforce_nor.out
 * rms_vel.out
 * viscf.out

Input variables in bou.in
-------------------------
```
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
bou.in
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
N1M             : Number of computational cells in the wall normal direction
N2M, N3M        : Number of computational cells in the homogenous/periodic directions
NSST            : if (NSST = 3) third order Runge Kutta integrator
                  if (NSST = 1) Adams Bashforth integrator
NREAD           : if (NREAD = 1) Start simulation by reading in continuation files
                  if (NREAD = 0) Start simulation without continuation files, initial
                                 conditions created in the code 

NTST            : Total number of time steps
WALLTIMEMAX     : Maximum wall time of the simulation
IRESET          : if (IRESET = 1) Time in simulation resets to 0
                  if (IRESET = 0) Time in simulation unchanged

XLEN, YLEN, ZLEN: Length of periodic directions 

REN             : Reynolds number
DT              : Time step to start the simulation
RESID           : Maximum allowable residual for mass conservation
CFLMAX          : Maximum CFL number for simulation (set to 1.2)

TSTA            : Time to start the statistics
STAREAD         : if (STAREAD = 1) read in statistics from continuation files
                  if (STAREAD = 0) don't read statistics

IDTV            : if (IDTV = 1) variable time stepping
                  if (IDTV = 0) fixed time stepping
DTMAX           : Maximum time step for variable time stepping

MOVIE           : if (MOVIE = 1) flow field is saved for visualization, h5 files are stored in
                  directory flowmov
                  if (MOVIE = 0) flow field is not saved
DUMPMOVIET      : Timestep to save the flow field

VORTR           : Vortex core radius
VORTSPEED:      : Vortex circulation

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
```

Input Variables in mlspart.in
-------------------------
```
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
mlspart.in
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
IMLSFORC        : if (IMLSFORC = 1) MLS forcing onto the flow is turned on
                  if (IMLSFORC = 0) MLS forcing onto the flow is turned off
IMLSSTRUC       : if (IMLSSTRUC = 1) MLS forcing onto the object is turned on
                  if (IMLSFORC = 0) MLS forcing onto the object is turned off
PREAD           : if (PREAD = 1) read in the position of the bodies rather than starting them from spos.in coordinates
                  if (PREAD = 0) position bodies at spos.in coordinates

SCLF            : Scale factor for the objects
RHOP            : Density of the particle

GTSFX           : GTS file/files to be read in
DATFX           : Name of restart file for the objects

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
```
