function [v_opt,flag,relres,iter,resvec] = SOA_GMRES_Hv_F(case_number, method, data_file)


% addpath(genpath('/Users/Ramsha/.../SWE Bath Data Assimilation 1D'))
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

G_eta = false;

if G_eta
    [F] = SOA_define_F(case_number, data_file);
else
    [F] = SOA_define_F_beta(case_number, data_file);
end
    

tol = 1e-8;
maxit = 256;
restart = 20;
% h = 6/256; X = (-3:h: 3-h)'; 
H = @(v) SOA_hess_bath(v, case_number, data_file); 
%v0 = 1e5*sin(10*X*pi/3) + 1e5*exp(- (X).^2);

if method == 1
    disp('commencing GMRES')
    [v_opt,flag,relres,iter,resvec] = gmres(H,F,restart,tol,maxit, [], [], []);
    
elseif method == 2
    [v_opt,flag,relres,iter,resvec] = bicgstab(H,F,[],maxit);
elseif method == 3
    disp('commencing BICGSTABl')
    [v_opt,flag,relres,iter,resvec] = bicgstabl(H,F,[],maxit);
end
    

