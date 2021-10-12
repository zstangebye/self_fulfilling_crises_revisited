# self_fulfilling_crises_revisited
Fortran Solution Code for "Self-Fulfilling Debt Crises Revisited"

There are three files written by the authors: globvar.f90, funcall.f90, and main.f90
The remaining files are open-source auxiliary files used in the solution

Code is written to be parallelized in MPI, but can be modified in globvar.f90 to allow for
serial execution.

Compile all auxiliary code first and then compile authors' code in order, e.g.,

mpif90  -fallow-argument-mismatch -ffree-line-length-none  blas.f cdfnor.f cumnor.f devlpl.f dinvnr.f ipmpar.f lbfgsb.f linpack.f mvndst.f spmpar.f stvaln.f timer.f globvar.f90 funcall.f90 main.f90

The flags "-fallow-argument-mismatch" and "-ffree-line-length-none" are technical and may not be necessary for all compilers

To run counterfactuals, simply alter the variables in the globvar file

1) To alter the likelihoods of belief regimes, adjust pes_bel and con_bel
2) To alter the share of auction revenue retained by the lenders in default, alter lend_share; when doing this one should also alter aminset, which governs the upper bound on the size of the grid, as it significantly increases the allowable indebtedness. We typically set aminset=-3.0 in this case
3) To alter the variance of the small shock, adjust epsilon_upper_bar
4) To accommodate long-term bonds, make the following changes, which both parameterize it correctly and facilitate convergence
  a) Increase sds to 0.07 (helps with convergence, see Chatterjee and Eyigungor [2012])
  b) Set lambdaset = .0638 
  c) Set aminset=-8.0
  d) Set con_bel = .0044
  e) Set defp0=.459 (this is the default cost)
  f) Increase A to 550 (Debt grid size)
  g) Set converge_tol to .001, as the iterations cycle below this tolerance threshold without affecting relevant moments
