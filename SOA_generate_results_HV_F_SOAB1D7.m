clear all


%%

 
for case_number = 1:3
    
if case_number == 1
    ampl = [0.01, 0.05, 0.07, 0.1, 0.15, 0.17, 0.2, 0.25,0.27, 0.3];
elseif case_number == 2
    ampl = [0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2];
elseif case_number == 3
    ampl = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1];
end
    
    
    
    for iter = 1:length(ampl)
        %% Set DA Data Path
        data_file = sprintf('Results/SOAB1D7/DA Results/SOAB1D7_case_%d_amplitude_%d.mat', case_number, iter);
        disp(sprintf('Case %d  iter %d commencing', case_number, iter));
        
        %% Run optimization for v s.t. Hv = F 
        [v_opt,flag,relres,iteration,resvec] = SOA_GMRES_Hv_F(case_number, 3, data_file);
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
        
        str =sprintf('Results/SOAB1D7/SOA Results/SOAB1D7_results_case_%d_amplitude_%d.mat', case_number, iter);
        save(str,'v_opt','eta_hat', 'u_hat', 'Hv', 'F', 'beta_optimum','x0_inds');
        disp(sprintf('Case %d  iter %d completing', case_number, iter));
        
        clear(data_file)
        clear('v_opt','eta_hat', 'u_hat', 'Hv', 'F', 'beta_optimum','x0_inds')
    end
end







