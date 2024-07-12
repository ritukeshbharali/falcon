
program = falcon

# subdirs = ../src/util  ../src/material ../src/model ../src/module ../src/femodel ../src/solver
subdirs = ../src/constraints \
          $(wildcard ../src/fem/*/) \
          $(wildcard ../src/io/*/) \
          ../src/materials \
          $(wildcard ../src/rom/*/) \
          $(wildcard ../src/solvers/*/) \
          ../src/steppers \
          ../src/util
          
include $(JIVEDIR)/makefiles/prog.mk
include $(JIVEDIR)/makefiles/packages/*.mk

# With external solvers Umfpack, Intel Pardiso, MUMPS

# Get root directories
MKLROOT         = /opt/intel/oneapi/mkl/2022.0.2
MUMPSROOT       = /home/ritukesh/libraries/MUMPS_5.5.0

# Libraries and their dependencies
SUITESPARSELIBS = umfpack amd suitesparseconfig
INTELMKLLIBS    = mkl_intel_lp64 mkl_intel_thread iomp5 mkl_core
MUMPSLIBS       = dmumps mumps_common 
DEPLIBS         = lapack openblas gfortran esmumps scotch scotcherr pord mpiseq pthread

# Include directories
SUITESPARSEINC  = /usr/include/suitesparse
INTELMKLINC     = $(MKLROOT)/include
MUMPSINC        = $(MUMPSROOT)/include
DEPINC          = $(MUMPSROOT)/PORD/include
MPIINC          = $(MUMPSROOT)/libseq

# Library directories
INTELMKLLIB     = $(MKLROOT)/lib/intel64
MUMPSLIB        = $(MUMPSROOT)/lib

MY_LIBS          = $(SUITESPARSELIBS) $(INTELMKLLIBS) $(MUMPSLIBS) $(DEPLIBS)
MY_INCDIRS       = ../src $(SUITESPARSEINC) $(INTELMKLINC) $(MUMPSINC) $(MPIINC) $(DEPINC) /usr/include
MY_LIBDIRS	     = . $(INTELMKLLIB) $(MUMPSLIB)  $(MUMPSROOT)/PORD/lib $(MPIINC) /usr/lib/
MY_CXX_STD_FLAGS = -DWITH_UMFPACK -DWITH_PARDISO -DWITH_MUMPS -fopenmp
MY_LD_FLAGS      = -static

# Without external solvers

#MY_LIBS          = 
#MY_INCDIRS       = ../src
#MY_LIBDIRS	     = 
#MY_CXX_STD_FLAGS = 
#MY_LD_FLAGS      = 
