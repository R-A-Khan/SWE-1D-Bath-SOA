function [v,fval,exitflag,output] = SOA_root_Hv_F(v_guess, case_number)


addpath(genpath('/Users/Ramsha/.../SWE Bath Data Assimilation 1D'))
addpath(genpath('/Users/Ramsha/... 2D'))


% Solves for optimal v minimising ||Hv-F||_inf in SOA implementation
% Input:
%       v_guess = initial guess for beta_hat (N x 1)
%
%       Results of optimal DA scheme
%       case_number = 1 ; Gaussian IC, Gaussian Bath   
%       case_number = 2 ; Gaussian IC, Sandbar Bath   
%       case_number = 3 ; Periodic IC, Gaussian Bath  
%
% Output:
%       v = local minimiser of (Hv-F, inf)


[F] = SOA_define_F(case_number);

options = optimset('Display','iter', 'Tolfun', 1e-8);
f = @(beta_hat) SOA_define_Hv_F(beta_hat, case_number, F);
v0 = v_guess;
[v, fval,exitflag,output] = fminunc(f,v0,options);
