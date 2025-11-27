External linear solvers
=======================

Wrapper classes have been written to solve the system of equations with external linear solvers. Available options are: 

- [AMGX (NVIDIA)](https://github.com/NVIDIA/AMGX)
- [Pardiso (Intel)](https://www.intel.com/content/www/us/en/develop/documentation/onemkl-developer-reference-c/top/sparse-solver-routines/onemkl-pardiso-parallel-direct-sparse-solver-iface.html).
- [MUMPS](https://mumps-solver.org/)
- [Pardiso (Panua)](https://panua.ch/pardiso/)
- [Umfpack (Suitesparse)](https://people.engr.tamu.edu/davis/suitesparse.html)

These libraries are covered by their own licenses, and the user is responsible for obtaining copies and installation. *Currently, only the sequential, multi-threaded version of falcon works with these libraries.*