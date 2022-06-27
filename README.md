# SWE-1D-Bath-SOA
A second order adjoint sensitivity analysis to quantify the sensitivity of wave forecasts derived using data assimilation, to placement and configuration of sensors. Full description of model derivation and evaluation can be found at:

Khan, R.A. and Kevlahan, N.K.-R., 2022. Variational Assimilation of Surface Wave Data for Bathymetry Reconstruction. Part II: Second Order Adjoint Sensitivity Analysis. Tellus A: Dynamic Meteorology and Oceanography, 74(1), pp.187–203. DOI: http://doi.org/10.16993/tellusa.36

## Usage
Adjust model parameters in 
* __SOA_generate_results.m__. 
 
Utilises the data assimilation algorithm and results obtained via thw SWE-1D-Bath project.

Verification of numerical solvers used in Second Order Adjoint scheme done using a kappa convergence test in 
* __Kappa_test_bath.m__

Visualize results using 
* __plots_hv_f.m__ 

Built using Matlab 2018a, including Optimization toolbox.

