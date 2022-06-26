% Data Assimilation Scheme for Bathymetry as Control Variable
 
clear all

sigma = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];

for case_number = 1
    for l = 1:length(sigma)
%%
% Fixed parameters for the data assimilation calculations
cfl         = 1/3;
n_obs       = 60;% Number of observation points
rand_x0     = false;        % Use random observation points
ntrial      = 1;            % Number of trials
iter_max    = 500  ;           % Max number of iterations
line_min    = true;        % Find optimal step
tau_n       = 8e3;          % Fixed step
nl          = true;         % nonlinear or linear equations
N           = 256 ;          % Number of grid points
xmax        = 3;            % domain size: xmin = -xmax
tmax        = 2*xmax;       % control time
smooth_grad = true;        % Smooth L2 gradient to Sobolev gradient
x0_min      =  -xmax;     % Location of first observation point
% dmu         = 0.05*xmax;     % Spacing between observation points
% dmu         = abs((xmax-1e-4) - x0_min)/n_obs;
 
%%%%%%%%% Parameters that we are varying %%%%%%
 
amp_bathy = 0.1;            % amplitude of bathymetry (0.3 of full depth is max for periodic ic)

if case_number == 1
    ic_per    = false;          % is initial condition periodic (or gaussian)?
    bathy_per = false;           % is bathymetry periodic (or gaussian)?
    bathy_np_type = 1;          % Gaussian = 1, Sandbar profile =2
    filt        = 0.05;          % filtering parameter for Sobolev gradient
elseif case_number == 2
    ic_per    = false;          % is initial condition periodic (or gaussian)?
    bathy_per = false;
    bathy_np_type = 2;          % Gaussian = 1, Sandbar profile =2
    filt        = 0.65;          % filtering parameter for Sobolev gradient
elseif case_number == 3
    ic_per    = true;          % is initial condition periodic (or gaussian)?
    bathy_per = false;
    bathy_np_type = 1;
    filt        = 0.03;          % filtering parameter for Sobolev gradient
end
 
amp_ic    = 0.01*amp_bathy; % amplitude of initial condition
k_ic      = 4;              % number of wavelenths of ic in domain
 
k_bathy   = 2;              % number of wavelengths of bathymetry in domain
phi_bathy = pi/2;         % phase shift of bathymetry with respect to initial condition
 
a         = 1;              % range of step size parameter, where the search for next step size is in [tau_n-tau_n/a,tau_n+tau_n/c]
c         =  1/10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
conj_grad_type = 0;
iter_chunk     = 10 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
noisy_obs  = false;  
print_iter = false;        % Print each iteration values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

disp(' ');
j_nl=1;

    nl=1;
    for k = 1: length(n_obs)
        dmu       = abs((xmax-1e-4) - x0_min)/n_obs(k);
         [eta_exct0(:,k), beta_optimum(:,k), beta_exct0(:,k), X, err(:,k), err_std(:,k), grad(:,k), grad_std(:,k), cost(:,k), cost_std, Y_ex(:,:,k), Y_opt(:,:,k), Y_adj_opt(:,:,k), x0_inds] = ...
           data_assimil_bath_2020(N,n_obs(k),rand_x0,x0_min,dmu, ntrial,tmax,xmax,iter_max,line_min,smooth_grad,filt, tau_n, amp_bathy, ic_per, bathy_per, ...
    bathy_np_type, amp_ic, k_ic, k_bathy, phi_bathy, sigma(l), conj_grad_type, iter_chunk, a, c, cfl, noisy_obs, print_iter);
    end
    disp(' ');

%% Saving Data
% 
% str =sprintf('SOAB1D5_case_%d_nObs_%d.mat', case_number, n_obs);
% 
% save(str, 'ntrial', 'n_obs', 'rand_x0','x0_min', 'dmu', 'line_min','smooth_grad', 'iter_max', ...
%     'eta_exct0', 'beta_optimum', 'beta_exct0', 'X', 'err', 'err_std', 'grad', 'grad_std', 'cost', 'cost_std', 'Y_ex', 'Y_opt', 'Y_adj_opt', 'x0_inds');
% 
% disp(sprintf('Case %d iteration %d with n_obs = %d completed', case_number, l, n_obs));
    end
end