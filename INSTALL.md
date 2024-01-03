# Installation

falcon is a Finite Element Analysis package based on the [Jem-Jive](https://software.dynaflow.com/jive/) libraries. As such, one needs to obtain and compile the Jem and Jive libraries. Prior to that, some pre-requisites must be installed.

## Pre-requisites
> (sudo) apt-get install build-essential g++ zlib1g-dev libreadline-dev freeglut3-dev

## Compile Jem-Jive libraries
  - Obtain the libraries from the [Dynaflow website](https://software.dynaflow.com/jive/)
  - Unpack archives Jem and Jive (say, in /home/user/libraries/jemjive)
  - Enter jem directory and './configure'
  - *(for MPI version)   ./configure --with-mpi --mpi-incdirs=/usr/include/x86_64-linux-gnu/openmpi --mpi-libdirs=/usr/lib/x86_64-linux-gnu/openmpi/lib --mpi-libs="mpi_cxx mpi"*
  - 'make' the library
  - export JEMDIR='/home/user/libraries/jemjive/jem-3.0'
  - Enter Jive directory and './configure', and 'make' the library
  - export JIVEDIR='/home/user/libraries/jemjive/jive-3.0'

## External libraries
falcon includes wrapper classes to solve the system of equations with external linear solvers. Available options are: 

- [AMGCL](https://github.com/ddemidov/amgcl)
- [AMGX (NVIDIA)](https://github.com/NVIDIA/AMGX)
- [Pardiso (Intel)](https://www.intel.com/content/www/us/en/develop/documentation/onemkl-developer-reference-c/top/sparse-solver-routines/onemkl-pardiso-parallel-direct-sparse-solver-iface.html).
- [MUMPS](https://mumps-solver.org/)
- [Pardiso (Panua)](https://panua.ch/pardiso/)
- [Umfpack (Suitesparse)](https://people.engr.tamu.edu/davis/suitesparse.html)

These libraries are covered by their own licenses, and the user is responsible for obtaining copies and installation. *Currently, only the sequential, multi-threaded version of falcon works with these libraries.*

## Build falcon
  - Obtain the source files
  - Enter the 'build' directory
  - Build optimized executable with 'make opt'
  - Enter the 'tests' directory, execute 'do_tests.sh' script
  - Larger (field-specific) problems may be run from 'demo' directory