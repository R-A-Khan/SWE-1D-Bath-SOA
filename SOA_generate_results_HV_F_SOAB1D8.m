clear all
cd(('/Users/Ramsha/Dropbox/McMaster University/PhD Research/Matlab Code/SOA_Implementation'))

%%

for case_number = 1
    

        %% Set DA Data Path
        data_file = sprintf('Results/SOAB1D3/DA Results/SOAB1D3_case_%d_x0_1.mat', case_number);
        disp(sprintf('Case %d commencing', case_number));
        
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
        
        str =sprintf('Results/SOAB1D8/SOA Results/SOAB1D8_results_case_%d_x0_1_GMRES.mat', case_number);
        save(str,'flag','relres','iteration','resvec','v_opt','eta_hat', 'u_hat',...
            'Hv', 'F', 'beta_optimum','x0_inds');
        disp(sprintf('Case %d completing', case_number));
        
        clear(data_file)
        clear('v_opt','eta_hat', 'u_hat', 'Hv', 'F', 'beta_optimum','x0_inds')
end






