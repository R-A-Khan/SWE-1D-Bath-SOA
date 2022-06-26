function [F] = SOA_define_F(case_number, data_file)
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

% if case_number == 1
%     if ~double_obs 
%     load('Bath/Case DAB1D21/Case_DAB1D21i.mat')
% %     disp('Single sided observation points');
%     else
%     load('Bath/Case DAB1D22/Case_DAB1D22i.mat')
% %      disp('Double sided observation points');
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

%%  Define G(eta, u, lambda) 

% G is in Y x Y_p 
% implicitly depends on x
eta_exct = Y_ex(1:N,:);     % (x,t)
u_exct = Y_ex(N+1:2*N,:);
G = trapz( T, (eta_opt - eta_exct).^2, 2) ;

%% Define d/dx for eta, u

% eta:= eta(x,t), d_eta/dt:= d_eta/dt(x,t) 
d_eta_dx = zeros(N,Nt);

    for col = 1:Nt
        % boundary conditions
        d_eta_dx(N,col) = ( eta_opt(1,col) - eta_opt(N-1,col) )/(2*delta_x);
        d_eta_dx(1,col) =  ( eta_opt(2,col) - eta_opt(N,col)   )/(2*delta_x);
        
        for row = 2:N-1
            d_eta_dx(row,col) = ( eta_opt(row+1,col) - eta_opt(row-1,col) )/(2*delta_x);
        end
    end
 
% G:= G(x), dG/dx:= dG/dx(x)
dG_dx = zeros(N,1);
dG_dx(N) = ( G(1) - G(N-1) )/(2*delta_x);
dG_dx(1) = ( G(2) - G(N)   )/(2*delta_x);

for row = 2:N-1
    dG_dx(row) = ( G(row+1) - G(row-1) )/(2*delta_x);
end


% u:= eta(x,t), du/dt:= du/dt(x,t) 
du_dx = zeros(N,Nt);

    for col = 1:Nt
        % boundary conditions
        du_dx(N,col) = ( u_opt(1,col) - u_opt(N-1,col) )/(2*delta_x);
        du_dx(1,col) =  ( u_opt(2,col) - u_opt(N,col)   )/(2*delta_x);
        
        for row = 2:N-1
            du_dx(row,col) = ( u_opt(row+1,col) - u_opt(row-1,col) )/(2*delta_x);
        end
    end
    

    
%% Define dG wrt eta and u   
d_eta_dx(:,1) = 0;
d_eta_du   = d_eta_dx ./( edge2mid_2D_x(du_dx) );
d_eta_du(isnan(d_eta_du)) = 0;



dG_deta = trapz(T, -2*( eta_exct - eta_opt )  , 2);
dG_du   = trapz(T, -2*d_eta_du.*( eta_exct - eta_opt ), 2);



%% Solve forced adjoint equation 

% Find gamma, psi
bw = @forced_FOA_bath;
H_foa = zeros(2*N,1);
[~, ~, ~, Y_foa] = BW_solve_forcedFOA(delta_x, dt, H_foa, 0, tmax,  bw, lambda, u_opt, eta_opt, dG_deta, dG_du);
Y_foa = fliplr(Y_foa);
gamma = Y_foa(1:N,:);
psi   = Y_foa(N+1:2*N,:);
% figure(1);surf(gamma, 'edgecolor','none')
% disp(['Inf norm of gamma = ', sprintf('%0.4f', norm(gamma, inf))])


%% Define dG/dlambda
% dg/dl(x,t) = dg/d_eta(x,t) * ( d_eta/dx(x,t)  ./ dl/dx(x) )
dl_dx = zeros(N,1);
dl_dx(N) = ( beta_optimum(1) - beta_optimum(N-1) )/(2*delta_x);
dl_dx(1) = ( beta_optimum(2) - beta_optimum(N)   )/(2*delta_x);

for row = 2:N-1
    dl_dx(row) = ( beta_optimum(row+1) - beta_optimum(row-1) )/(2*delta_x);
end

% if adjust
% dl_dx(1:150) = 0;
% dl_dx(225:end) = 0;
% end
% 


dG_dl = dG_dx./dl_dx;

%% Define F
d_gamma_dx= zeros(N,Nt);
for j = 1:Nt
    d_gamma_dx(:,j) = cent_diff_u(delta_x,gamma(:,j) );
end

prod = u_opt.*d_gamma_dx;
int_F = trapz(T,prod,2);

F = dG_dl - int_F;

end