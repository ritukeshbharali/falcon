External linear solvers
=======================
The following wrapper classes have been written to accommodate external
solvers for the linear system of equations:

- [AMGCL](https://amgcl.readthedocs.io/en/latest/)
- [Intel Pardiso](https://www.intel.com/content/www/us/en/develop/documentation/onemkl-developer-reference-c/top/sparse-solver-routines/onemkl-pardiso-parallel-direct-sparse-solver-iface.html).
- [MUMPS](https://mumps-solver.org/)
- [Umfpack](https://people.engr.tamu.edu/davis/suitesparse.html)

These libraries are covered by their own licenses, and the user is responsible for obtaining copies and installation. *Currently, only the sequential, multi-threaded version of falcon works with these libraries.*
