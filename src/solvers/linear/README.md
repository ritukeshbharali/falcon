External solvers
================
The following wrapper classes have been written to accomodate external
solvers for the linear system of equations:

- [Umfpack](https://people.engr.tamu.edu/davis/suitesparse.html)
- [Intel Pardiso](https://www.intel.com/content/www/us/en/develop/documentation/onemkl-developer-reference-c/top/sparse-solver-routines/onemkl-pardiso-parallel-direct-sparse-solver-iface.html).
- [MUMPS](http://mumps.enseeiht.fr/)

These libraries are covered by their own licenses, and the user is responsible for obtaining copies and installation. *Currently, only the sequential, multi-threaded version of jiveFEA works with these libraries.*
