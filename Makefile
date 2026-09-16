# ===========================================================================
#  Makefile for VP_PIC
#
#    make            build exe/VP_PIC with optimization
#    make debug      rebuild with run-time checks and backtraces
#    make run        run it (see ARGS below)
#    make clean      remove the objects and the executable
#    make veryclean  also remove exe/ with everything in it
#    make help       this list
#
#  Variables that can be set on the command line, e.g. "make FC=ifort":
#
#    FC           Fortran compiler (gfortran or ifort).
#    HDF5_INC     -I flag for the HDF5 Fortran module.
#    HDF5_LIBS    -L and -l flags for the HDF5 libraries.
#    OMP_THREADS  Threads used by "make run".
#    ARGS         Arguments passed by "make run", for instance
#                 make run ARGS="base.par Nt=2000 directory=run7"
# ===========================================================================

SRCDIR := src
OBJDIR := objs
EXEDIR := exe
EXE    := $(EXEDIR)/VP_PIC


# ---------------------------------------------------------------------------
#  Compiler and flags
# ---------------------------------------------------------------------------

# "make" itself predefines FC, so only replace it when the user has not.

ifeq ($(origin FC),default)
  FC := gfortran
endif

# Flags are split so that the optimized and the debug build share everything
# except the optimization level and the checks.
#
#   -ffree-form / -free            free form source
#   -fopenmp / -qopenmp            the parallel loops
#   -J / -module                   where to leave the .mod files
#   -fallow-argument-mismatch      the HDF5 Fortran interface passes buffers
#                                  of different types to the same argument
#   -Wno-unused-dummy-argument     the B-spline shape functions share one
#                                  signature, so some ignore an argument

ifeq ($(FC),gfortran)
  FFLAGS_BASE  := -ffree-form -fopenmp -fallow-argument-mismatch -J$(OBJDIR) \
                  -Wall -Wno-unused-dummy-argument
  FFLAGS_OPT   := -O3 -funroll-loops
  FFLAGS_DEBUG := -O0 -g -fcheck=all -fbacktrace
else ifeq ($(FC),ifort)
  FFLAGS_BASE  := -free -shared-intel -qopenmp -warn -module $(OBJDIR)
  FFLAGS_OPT   := -O3 -xhost -align array64byte
  FFLAGS_DEBUG := -O0 -g -check all -traceback -fpe0
else
  $(error Unknown Fortran compiler "$(FC)". Add its flags to the Makefile.)
endif

# BUILD is set to "debug" by the debug target below. FFLAGS is expanded when
# used, not here, so that the target can still change BUILD.
#
# The debug build checks bounds and prints backtraces, but does NOT trap
# floating point exceptions: initial_data builds the action-angle variables
# for every cell of the (r,p) box, and the cells that are unbound produce a
# NaN on purpose, which is what marks the particle for removal.

BUILD  ?= opt
FFLAGS  = $(FFLAGS_BASE) $(if $(filter debug,$(BUILD)),$(FFLAGS_DEBUG),$(FFLAGS_OPT))

# HDF5 (Fortran bindings, serial build). These are the paths of a Debian or
# Ubuntu libhdf5-dev; "h5fc -show" prints the ones of any other install. The
# rpath records the library directory in the executable, so it runs without
# LD_LIBRARY_PATH.

HDF5_INC  ?= -I/usr/include/hdf5/serial
HDF5_LIBS ?= -L/usr/lib/x86_64-linux-gnu/hdf5/serial \
             -Wl,-rpath,/usr/lib/x86_64-linux-gnu/hdf5/serial \
             -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5

# Threads for "make run". OMP_PLACES/OMP_PROC_BIND put one thread per physical
# core, instead of two sharing the hyperthreads of the same core.

OMP_THREADS ?= 4
ARGS        ?=


# ---------------------------------------------------------------------------
#  Sources
# ---------------------------------------------------------------------------

# The modules go first: everything else uses them. Listed by hand because the
# order matters, the rest is picked up automatically.

MODULES := parameters paramfile arrays utils functions hdf5_io raw_io

MODOBJS := $(addprefix $(OBJDIR)/,$(addsuffix .o,$(MODULES)))
ALLOBJS := $(patsubst $(SRCDIR)/%.f90,$(OBJDIR)/%.o,$(wildcard $(SRCDIR)/*.f90))
OBJS    := $(filter-out $(MODOBJS),$(ALLOBJS))

# The first rule in a makefile is the one "make" builds by default, and the
# dependencies just below are rules, so the default goal is stated explicitly.

.DEFAULT_GOAL := all

# Which module uses which, so that "make -j" still compiles them in an order
# that works.

$(OBJDIR)/paramfile.o : $(OBJDIR)/parameters.o
$(OBJDIR)/utils.o     : $(OBJDIR)/parameters.o $(OBJDIR)/arrays.o
$(OBJDIR)/hdf5_io.o   : $(OBJDIR)/parameters.o $(OBJDIR)/arrays.o
$(OBJDIR)/raw_io.o    : $(OBJDIR)/parameters.o $(OBJDIR)/arrays.o
$(OBJS)               : $(MODOBJS)


# ---------------------------------------------------------------------------
#  Rules
# ---------------------------------------------------------------------------

.PHONY: all debug run clean veryclean help

all: $(EXE)

$(OBJDIR)/%.o : $(SRCDIR)/%.f90 | $(OBJDIR)
	@ echo "COMPILING FILE: $(notdir $<)"
	@ $(FC) $(FFLAGS) $(HDF5_INC) -I$(OBJDIR) -c $< -o $@

$(EXE) : $(MODOBJS) $(OBJS) | $(EXEDIR)
	@ echo
	@ echo "LINKING ..."
	@ $(FC) $(FFLAGS) $^ $(HDF5_LIBS) -o $@
	@ echo
	@ echo "COMPILATION DONE!"
	@ echo

$(OBJDIR) $(EXEDIR) :
	@ mkdir -p $@

# A debug build cannot reuse optimized objects, so start from scratch.

debug :
	@ $(MAKE) --no-print-directory clean
	@ $(MAKE) --no-print-directory BUILD=debug

run : $(EXE)
	@ cd $(EXEDIR) && OMP_NUM_THREADS=$(OMP_THREADS) OMP_PLACES=cores \
	  OMP_PROC_BIND=close ./VP_PIC $(ARGS)

clean :
	@ /bin/rm -rf $(OBJDIR) $(EXE)

# Deletes exe/ whole, with every run that was left inside it.

veryclean :
	@ /bin/rm -rf $(OBJDIR) $(EXEDIR)

help :
	@ echo
	@ echo "make            Build $(EXE) with optimization"
	@ echo "make debug      Rebuild with run-time checks and backtraces"
	@ echo "make run        Run it with OMP_NUM_THREADS=$(OMP_THREADS)"
	@ echo "                  make run ARGS=\"base.par Nt=2000 directory=run7\""
	@ echo "                  make run OMP_THREADS=8"
	@ echo "make clean      Delete $(OBJDIR)/ and $(EXE)"
	@ echo "make veryclean  Delete $(OBJDIR)/ and $(EXEDIR)/ with ALL its contents"
	@ echo
	@ echo "Compiler: FC=$(FC)"
	@ echo
