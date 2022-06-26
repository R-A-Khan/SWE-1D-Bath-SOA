function [A_inf] = SOA_define_Hv_F(beta_hat, case_number, F)

addpath(genpath('/Users/Ramsha/Dropbox/McMaster University/PhD Research/Matlab Code/SWE Bath Data Assimilation 1D'))
addpath(genpath('/Users/Ramsha/Dropbox/McMaster University/PhD Research/Matlab Code 2D'))

%% Find Hessian of J given beta_hat

[~, Hess_mesh] = SOA_hess_bath( beta_hat, case_number);


%% Define A
A = Hess_mesh - F;
A_inf = norm(A, inf);
% figure;surf(A, 'edgecolor','none')
% disp(['Inf norm of A = ', sprintf('%0.4f', A_inf)])

