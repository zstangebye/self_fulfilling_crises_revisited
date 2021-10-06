# self_fulfilling_crises_revisited
Fortran Solution Code for "Self-Fulfilling Debt Crises Revisited"

There are three files written by the authors: globvar.f90, funcall.f90, and main.f90
The remaining files are open-source auxiliary files used in the solution

Code is written to be parallelized in MPI, but can be modified in globvar.f90 to allow for
serial execution.

Compile all auxiliary code first and then compile authors code in order: globvar.f90 funcall.f90 main.f90
