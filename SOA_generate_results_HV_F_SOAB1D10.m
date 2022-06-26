clear all
%cd(('/Users/Ramsha/Dropbox/McMaster University/PhD Research/Matlab Code/SOA_Implementation'))

%%

for case_number = 1:3
          
for l = 1:4
        disp(sprintf('Case %d iteration %d commencing', case_number, l));  
        %% Set DA Data Path
        data_file = sprintf('Results/SOAB1D10/DA Results/SOAB1D10_case_%d_x0_%d.mat', case_number,l);

        %% Run optimization for v s.t. Hv = F 
        [v_opt,flag,relres,iteration,resvec] = SOA_GMRES_Hv_F(case_number, 1, data_file);
        [Hv] = SOA_hess_bath( v_opt, case_number, data_file);
        [F] = SOA_define_F(case_number, data_file);
        load(data_file)
        
    
        %% Get DG/dm
        N = length(beta_exct0);
        eta_opt = Y_opt(1:N,:);     % (x,t)
        u_opt   = Y_opt(N+1:2*N,:); % (x,t)
        lambda = beta_optimum;      % (x)

        cfl = 1/3;
        xmax = 3; 
        xmin = -xmax;
        tmax = 6;
        tmin = 0;
        delta_x = (xmax - xmin)/N;
        dt = delta_x*cfl;
        T = tmin:dt:tmax;
        fw_tangent_BATH = @fw_SWE_NL_tangent_model_BATH;


        u = zeros(N,1);
        eta_hat0 = zeros(N,1);
        H_hat = [eta_hat0 ; u];
        [~, ~, eta_hat_x0_pts, ~, T_hat, Y_hat] = FW_solve_tangent_BATH(delta_x, dt, H_hat, 0, tmax, u_opt, eta_opt, beta_optimum, v_opt, x0_inds, fw_tangent_BATH);
        eta_hat = Y_hat(1:N,:);
        u_hat = Y_hat(N+1:2*N,:);
 

        %% Save Results
        
        str =sprintf('Results/SOAB1D10/SOA Results/SOAB1D10_results_case_%d_N_obs_%d_GMRES.mat', case_number,l);
        save(str,'flag','relres','iteration','resvec','v_opt','eta_hat', 'u_hat',...
            'Hv', 'F', 'beta_optimum','x0_inds');
        
        
        clear(data_file)
        clear('v_opt','eta_hat', 'u_hat', 'Hv', 'F', 'beta_optimum','x0_inds')
end
disp(sprintf('Case %d completing', case_number));
end





