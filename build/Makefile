#---------------------------------------------------------
#
# General makefile for Nbody_acc; based on Phantom
#
#---------------------------------------------------------

.KEEP_STATE:

KNOWN_SYSTEM=no
SHELL=/bin/bash
VPATH="${RUNDIR}" ../src/main ../src/setup 
BINDIR= ../bin

OPENMP=yes
DOUBLEPRECISION=yes
ENDIAN=little 
DEBUG=no 
LDFLAGS=
FPPFLAGS=

ifeq ($(SYSTEM), gfortran)
    FC= gfortran
    FFLAGS+= -O3 -Wall -Wno-unused-dummy-argument -frecord-marker=4 -gdwarf-2 \
             -finline-functions-called-once -finline-limit=1500 -funroll-loops -ftree-vectorize \
             -std=f2008 -fall-intrinsics
    DBLFLAG= -fdefault-real-8 -fdefault-double-8
    DEBUGFLAG= -g -fcheck=all -ffpe-trap=invalid,zero,overflow -finit-real=nan -finit-integer=nan -fbacktrace
    ENDIANFLAGBIG= -fconvert=big-endian
    ENDIANFLAGLITTLE= -fconvert=little-endian
    OMPFLAGS = -fopenmp
    KNOWN_SYSTEM= yes
endif

ifeq ($(SYSTEM), ifort)
    FC= ifort
    FFLAGS= -O3 -inline-factor=500 -mcmodel=medium -shared-intel -warn uninitialized -warn unused -warn truncated_source -no-wrap-margin
    DBLFLAG= -r8
    DEBUGFLAG= -check all -WB -traceback -g -debug all # -fpe0 -fp-stack-check -debug all -noarg_temp_created
    ENDIANFLAGBIG= -convert big_endian
    ENDIANFLAGLITTLE= -convert little_endian
    OMPFLAGS= -qopenmp
    KNOWN_SYSTEM= yes 
endif

ifneq ($(KNOWN_SYSTEM), yes)
	@echo "Please set SYSTEM variable in your bashrc"
	@${MAKE} err;
endif

#
# Enter the setup file 
#
SETUPFILE= setup_turbcore.f90

#
# Configuration options - physics 
#
BINARY= no
GRAVITY= yes

ifeq ($(GRAVITY), yes)
    FPPFLAGS += -DGRAVITY
endif
ifeq ($(BINARY), yes)
    FPPFLAGS += -DBINARY
endif

#
# Configuration options - software 
#
ifeq ($(OPENMP), yes)
    FFLAGS += ${OMPFLAGS}
endif
ifeq ($(DOUBLEPRECISION), yes)
    FFLAGS += ${DBLFLAG}
endif
ifeq ($(DEBUG), yes)
    FFLAGS += ${DEBUGFLAG}
endif 
ifeq ($(ENDIAN), big)
    FFLAGS += ${ENDIANFLAGBIG}
endif
ifeq ($(ENDIAN), little)
    FFLAGS += ${ENDIANFLAGLITTLE}
endif


#
# define the implicit rule to make a .o file from a .f90 file
#
.SUFFIXES:
.SUFFIXES: .o .f90 .F90
%.o : %.f90
	$(FC) -c $(FFLAGS) $< -o $@
%.o : %.F90
	$(FC) -c $(FFLAGS) ${FPP_PREFIX} $(FPPFLAGS) $< -o $@

#
# Main 
#
SOURCES= physcon.f90 units.f90 accrete.F90 timestep.f90 energies.f90 forces.F90 step_RK4.F90 \
         ptmass.F90 readwrite_dumpfile.F90 readwrite_infile.f90 evolve.F90 initial.F90  

OBJECTS1 = $(SOURCES:.f90=.o)
OBJECTS = $(OBJECTS1:.F90=.o)

.PHONY: nbodyacc

nbodyacc: $(OBJECTS) nbodyacc.o
	$(FC) $(FFLAGS) $(FPPFLAGS) -o $(BINDIR)/$@ $(OBJECTS) nbodyacc.o $(LDFLAGS)

#
# Setup 
#
SRCSETUP= physcon.f90 units.f90 timestep.f90 energies.f90 ptmass.F90 forces.F90 readwrite_dumpfile.f90 \
         vel_turbcube.f90 readwrite_infile.F90 ${SETUPFILE}   # ADD ALL FILES NEEDED BY SETUP 

OBJSETUP1= $(SRCSETUP:.f90=.o)
OBJSETUP= $(OBJSETUP1:.F90=.o) 

.PHONY: nbodyaccsetup

nbodyaccsetup: $(OBJSETUP) nbodyaccsetup.o
	$(FC) $(FFLAGS) $(FPPFLAGS) -o $(BINDIR)/$@  $(OBJSETUP) nbodyaccsetup.o $(LDFLAGS)

setup: nbodyaccsetup

#
# Clean 
#
clean: clean
	rm *.o; cd $(BINDIR); rm nbodyacc*; cd ../src/main; rm -f *.mod; cd ../setup; rm -f *.mod

err:
	$(error aborting);
