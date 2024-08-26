
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

# With external solver SuperLU_MT with pthreads (macOS)

# Get root directory
SUPERLUMTROOT    = /Users/ritukeshbharali/mylibs/superlu_mt-4.0.0/

# Libraries and their dependencies
SUPERLUMTLIBS    = superlu_mt_PTHREAD
DEPLIBS          = lapack openblas pthread

# Include directories
SUPERLUMTINC     = $(SUPERLUMTROOT)/SRC

# Library directories
SUPERLUMTLIB     = $(SUPERLUMTROOT)/lib

MY_LIBS          = $(SUPERLUMTLIBS) $(DEPLIBS)
MY_INCDIRS       = ../src $(SUPERLUMTINC) /usr/include
MY_LIBDIRS	     = . $(SUPERLUMTLIB) /usr/lib/
MY_CXX_STD_FLAGS = -DWITH_SUPERLUMT -DOSX
MY_LD_FLAGS      = -static

# Without external solvers

#MY_LIBS          = 
#MY_INCDIRS       = ../src
#MY_LIBDIRS	      = 
#MY_CXX_STD_FLAGS = 
#MY_LD_FLAGS      = 
