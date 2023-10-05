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
 * update_halo 2x
 * CalculateLocalDivergence
 * SolvePressureCorrection
 * update_halo
 * CorrectVelocity
 * CorrectPressure
 * update_halo 5x
 * CalcMLSForce
 * MLSForceVelocity

List of input and output files
------------------------------

Input files (.in)
 *  bou.in - Detailed information on the input variables is given in the next section
 *  mlspart.in - Input variables for the immersed bodies
 *  spos.in - Initial position of the body in the computational domain
 *  ppart.in
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

