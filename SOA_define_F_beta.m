function [F] = SOA_define_F_beta(case_number, data_file)
addpath(genpath('/Users/Ramsha/Dropbox/McMaster University/PhD Research/Matlab Code'))
% Defines the RHS function F in Hv = F for SOA implementation
% Input:
%       Results of optimal DA scheme
%       case_number = 1 ; Gaussian IC, Gaussian Bath   
%       case_number = 2 ; Gaussian IC, Sandbar Bath   
%       case_number = 3 ; Periodic IC, Gaussian Bath  
%
% Output:
%       F = dG/dl - int_0^T [ u(x,t) * dGamma/dx(x,t)] dt


%% Retrieve dataset for optimal DA results

load(data_file);

%% Extract Variables from Optimal DA Results

N = length(beta_exct0);
eta_opt = Y_opt(1:N,:);     % (x,t)
u_opt   = Y_opt(N+1:2*N,:); % (x,t)
lambda = beta_optimum;      % (x)

% get size of T
[~, Nt] = size(eta_opt);
cfl = 1/3;
xmax = 3; 
xmin = -xmax;
tmax = 6;
tmin = 0;
delta_x = (xmax - xmin)/N;
dt = delta_x*cfl;
T = tmin:dt:tmax;
adjust = false;

%%  Define G( lambda) 

% G is in Y x Y_p 
% implicitly depends on x
eta_exct = Y_ex(1:N,:);     % (x,t)
u_exct = Y_ex(N+1:2*N,:);
G =  (beta_exct0 - beta_optimum).^2;



%% Define dG/dlambda
dG_dl = zeros(N,1);
dG_dl = -2*(beta_exct0 - beta_optimum);




%% Define F
F = dG_dl ;

