
function [Hv] = SOA_hess_bath( v, case_number,data_file)

% addpath(genpath('/Users/Ramsha/.../SWE Bath Data Assimilation 1D'))
addpath(genpath('/Users/Ramsha/... 2D'))

% if case_number == 1
%     if ~double_obs 
%     load('Bath/Case DAB1D21/Case_DAB1D21i.mat')
%     else
%     load('Bath/Case DAB1D22/Case_DAB1D22i.mat')
%     end
% 
% elseif case_number == 2
%     if ~double_obs
%         load('Bath/Case DAB1D21/Case_DAB1D21ii.mat')
%     else
%         load('Bath/Case DAB1D22/Case_DAB1D22ii.mat')
%     end
% elseif case_number == 3
%     if ~double_obs
%         load('Bath/Case DAB1D21/Case_DAB1D21iii.mat')
%     else
%         load('Bath/Case DAB1D22/Case_DAB1D22iii.mat')
%     end
% else
% %     disp('You have entered an invalid case number')
% end

load(data_file)

%%

fw_tangent_BATH = @fw_SWE_NL_tangent_model_BATH;
bw_SOA_BATH = @bw_SWE_SOA_BATH;


% Extract variables from optimality system
N = length(beta_exct0);
eta_opt = Y_opt(1:N,:);     % (x,t)
u_opt   = Y_opt(N+1:2*N,:); % (x,t)
[~, Nt] = size(eta_opt);
cfl = 1/3;
xmax = 3; 
xmin = -xmax;
tmax = 6;
tmin = 0;
delta_x = (xmax - xmin)/N;
dt = delta_x*cfl;
T = tmin:dt:tmax;
nl = true;

% Optimal Adjoint Variables
% KEEP THESE IN TERMS OF TAU; that is how they are input in BW_SOA
eta_adj_opt = Y_adj_opt(1:N,:);
u_adj_opt   = Y_adj_opt(N+1:2*N,:);


    
      
% Solve tangent model for perturbed system beta + ep*beta_hat
u = zeros(N,1);
eta_hat0 = zeros(N,1);
H_hat = [eta_hat0 ; u];
[~, ~, eta_hat_x0_pts, ~, T_hat, Y_hat] = FW_solve_tangent_BATH(delta_x, dt, H_hat, 0, tmax, u_opt, eta_opt, beta_optimum, v, x0_inds, fw_tangent_BATH);
eta_hat = Y_hat(1:N,:);
u_hat = Y_hat(N+1:2*N,:);

% Solve SOA model for u_bar
H_bar = zeros(2*N,1);
[~, T_bar, Y_bar] = BW_solve_SOA_BATH(delta_x, dt, H_bar, 0, tmax, x0_inds, bw_SOA_BATH, eta_hat_x0_pts, u_opt, eta_opt, u_adj_opt, eta_adj_opt, u_hat, eta_hat, beta_optimum, v, T_hat, T, nl);


% Define Hessian

Y_adj_flip = fliplr(Y_adj_opt);
d_eta_adj_x_t= zeros(N,length(T));
for j = 1:length(T)
    d_eta_adj_x_t(:,j) = cent_diff_u(delta_x,(Y_adj_flip(1:N,j)));
end

Y_bar = fliplr(Y_bar);
d_eta_bar_x_t= zeros(N,length(T_bar));
for j = 1:length(T)
    d_eta_bar_x_t(:,j) = cent_diff_u(delta_x,(Y_bar(1:N,j)));
end
prod1 = u_hat.*d_eta_adj_x_t;
prod2 = u_opt.*d_eta_bar_x_t;
Hv = trapz(T,(prod1+prod2),2);


% figure(1); plot(Hessian_J);
% figure(2); surf(Hess_mesh, 'edgecolor', 'none')
% figure(3); plot(Int_Hess)
%     

