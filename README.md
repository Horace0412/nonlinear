# nonlinear

The problem that I am solving is a combinatorial programming problem. Here’s the mathematical structure: 
The choice variable is m=[m1, m2, m3, …, mN, m1m2, m1m3, …, m1mN, m2m3, m2m4, …, m3m4, …, mN-1mN], i.e.,  .
Essentially, only m1, …, mN are variables to choose. 

The problem takes as given two matrices, mu153×1 and sigma153×153, and compute two intermediate variables, mux = m @ mu and sigma2x = m @ sigma @ m’. sigma is a covariance matrix so it is positive semi-definite. sigma2x is computed using Cholesky decomposition as shown in sigma2_exp_gen.

The objective function is a nonlinear function of mux and sigma2x. I use a tabular method: I discretize mux and sigma2x into grid points, muy_grid and sigma2y_grid, and compute objective function value on grid points, profits_grid2. Then, in function gurobi_muSigma_sample, for a given m, I compute mux and sigma2x, then I do bilinear interpolation (i.e., using SOS2 constraints) to obtain objective function value. 
