#Makefile for parallel compiling
#=====================================================
# Note: use -C=all to catch array out of bounds errors
#============================================================================ 
# Compiler options
#============================================================================
# IFORT
#FC = h5pfc -r8 -O3 -fpp -ip -ipo -xAVX -axCORE-AVX2

#FC += -openmp
#FC += -g -traceback -check bounds
#FC += -debug all -warn all -check all -g -traceback
#FC += -fpe0 -ftrapuv


# GNU
#FC = h5pfc -O0 -g -cpp -Wextra -fdefault-real-8 -fbounds-check -fdefault-double-8 -fbacktrace
FC = h5pfc -r8 -ip -ipo -O3 -fpp -xHost -qno-openmp
#FC = mpif90 -O3 -cpp -Wunused-variable -Wunused-dummy-argument -fdefault-real-8 -fbounds-check -fdefault-double-8 -I/usr/include/hdf5/openmpi
#FC += -fopenmp

# IBM
#FC = mpixlf77_r -O3 -qautodbl=dbl4 -WF,-qfpp
#FC += -g -C -qfullpath -qinfo=all
#FC += -WF,-DMPI -WF,-DFREESLIP
#FC += -WF,-DSTATS
#FC += -WF,-DBALANCE
#FC += -WF,-DBLSNAP
#FC += -WF,-DDEBUG
#FC += -qstrict=all

# CRAY

# GENERAL
#FC += -DDEBUG
#FC += -DPLANEMOV
#FC += -DBLSNAP
#=======================================================================
# Library
#======================================================================

LINKSAFTER = -lfftw3  -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lhdf5_fortran -lhdf5  -lz -ldl -lm

#============================================================================ 
# make PROGRAM   
#============================================================================

PROGRAM = boutnp

OBJECTS = MLSForceVelocity.o findindices.o partindicesMLS.o tri_geo.o CalcMLSForce.o AuxiliaryRoutines.o \
          TriAuxRoutines.o MLSAuxRoutines.o mkmov_flo.o mkmov_str.o mlsStruc_closed_rigid.o write_3dfield.o\
          CalcMaxCFL.o CalcLocalDivergence.o CheckDivergence.o ExplicitTermsVX.o ExplicitTermsVY.o ExplicitTermsVZ.o \
          CreateGrid.o ReadInitCond.o InitArrays.o CreateInitialConditions.o ImplicitAndUpdateVX.o ImplicitAndUpdateVY.o ImplicitAndUpdateVZ.o  \
          OpenLogs.o CalcHITRandomForce.o matrix_transpose.o mpi_routines.o PointPartAuxRoutines.o \
          papero.o param.o CorrectPressure.o SolveTridX.o SolveTridY.o SolveTridZ.o  \
          StatRoutines.o allocate_mls_local.o allocate_trigeo.o HdfAuxRoutines.o QuitRoutine.o \
          PerTridiagSolve.o TimeMarcher.o CorrectVelocity.o CheckMaxVel.o InitPressSolv.o ReadMLSInput.o \
          SolvePressureCorrection.o CalcDissipation.o interp.o LocateLargeDivergence.o WriteRandForcCoef.o InitRandomForce.o CalcWriteQ.o \
	  write_to_vtk.o CalcVorticity.o MpiAuxRoutines.o InitTimeMarchScheme.o ReadInputFile.o \
	  InitPointP.o partvel_explicit.o UpdatePointParticlePosition.o CalcMaterialDerivative.o ReadPpartInput.o PointPartIO.o DumpMovieFiles.o

MODULES = param.o
#============================================================================ 
# Linking    
#============================================================================

$(PROGRAM) : $(MODULES) $(OBJECTS)
	$(FC) $(OBJECTS) $(LINKS)   -o $@  $(LINKSAFTER)

#============================================================================
#  Dependencies
#============================================================================

param.o: param.F90
	$(FC) -c $(OP_COMP) param.F90

calcforc.o:   calcforc.F90 $(MODULES)
	$(FC) -c $(OP_COMP) $< 

initfrc.o:   initfrc.F90 $(MODULES)
	$(FC) -c $(OP_COMP) $< 

%.o:   %.F $(MODULES)
	$(FC) -c $(OP_COMP) $< 

%.o:   %.F90 $(MODULES)
	$(FC) -c $(OP_COMP) $< 

um_ib.o:   um_ib.f90 $(MODULES)
	$(FC) -c $(OP_COMP) $< 
        

#============================================================================
#  Clean up
#============================================================================

clean :
	rm *.o *.mod boutnp

