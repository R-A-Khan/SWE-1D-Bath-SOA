for case_num = 2
test = [1 10 20];

    for k = 1:length(test)
str1 =sprintf('Results/SOAB1D3/DA Results/SOAB1D3_case_%d_x0_%d.mat', case_num, test(k));
str2 =sprintf('Results/SOAB1D3/SOA Results/SOAB1D3_results_case_%d_x0_%d.mat', case_num, test(k));

load(str1)
load(str2)

N = length(beta_exct0);
cfl = 1/3;
xmax = 3; 
xmin = -xmax;
tmax = 6;
tmin = 0;
delta_x = (xmax - xmin)/N;
dt = delta_x*cfl;
T = tmin:dt:tmax;


sens_at_obs = eta_hat(x0_inds,:);
sens_int_T = trapz(T, sens_at_obs, 2);
figure(k); 
semilogy(1:45,abs(sens_int_T), 'linewidth',2)
set(gca,'FontSize',14);
ylim([1e-7 1e-2])

xlabel('Observation Point','interpreter','latex')
ylabel('$ |\int_0^T \frac {\partial G}{\partial m} dt | $','interpreter','latex', 'fontsize',20)

str3 = sprintf('Results/SOAB1D3/Plots/case_%d_x0_%d_dG_dm_integrated.eps', case_num, test(k));
print(str3, '-depsc2')
str5 = sprintf('Results/SOAB1D3/Plots/case_%d_x0_%d_dG_dm_integrated.fig', case_num, test(k));
savefig(str5)

    end
end


%%
%% Plots for dG/dm integrated over time (fig 5 in SOA Tellus_
clear variables
for case_num = 1

    
 test =[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
    for k = [2 5 8]
str1 =sprintf('Results/SOAB1D6/DA Results/SOAB1D6_case_%d_sigma_%d.mat', case_num, k);
str2 =sprintf('Results/SOAB1D6/SOA Results/SOAB1D6_results_case_%d_sigma_%d.mat', case_num, k);
load(str1)
load(str2)

N = length(beta_exct0);
cfl = 1/3;
xmax = 3; 
xmin = -xmax;
tmax = 6;
tmin = 0;
delta_x = (xmax - xmin)/N;
dt = delta_x*cfl;
T = tmin:dt:tmax;


sens_at_obs = eta_hat(x0_inds,:);
sens_int_T = trapz(T, sens_at_obs, 2);

figure(k); 
semilogy(abs(sens_int_T), 'linewidth',2)
set(gca,'FontSize',14);



xlabel('Observation Point','interpreter','latex')
ylabel('$ |\int_0^T \frac {\partial G}{\partial m} dt |$','interpreter','latex', 'fontsize',20)

 ylim([1e-10 1e-3])

str3 = sprintf('Results/SOAB1D6/Plots/case_%d_sigma_%d_dG_dm_integrated.eps', case_num, k);
print(str3, '-depsc2')
str5 = sprintf('Results/SOAB1D6/Plots/case_%d_sigma_%d_dG_dm_integrated.fig', case_num, k);
savefig(str5)

    end
end


%%
%% Plots for dG/dm integrated over time
clear variables
for case_num = 1

    
 test =[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
    for k = [2 5 8]
str1 =sprintf('Results/SOAB1D5/DA Results/SOAB1D5_case_%d_sigma_%d.mat', case_num, k);
str2 =sprintf('Results/SOAB1D5/SOA Results/SOAB1D5_results_case_%d_sigma_%d.mat', case_num, k);
load(str1)
load(str2)

N = length(beta_exct0);
cfl = 1/3;
xmax = 3; 
xmin = -xmax;
tmax = 6;
tmin = 0;
delta_x = (xmax - xmin)/N;
dt = delta_x*cfl;
T = tmin:dt:tmax;


sens_at_obs = eta_hat(x0_inds,:);
sens_int_T = trapz(T, sens_at_obs, 2);

figure(k); 
semilogy(abs(sens_int_T), 'linewidth',2)
set(gca,'FontSize',14);
str4 = sprintf('Bathymetry Gaussian with Standard Deviation %0.2f',test(k));
title(str4, 'interpreter', 'latex')


xlabel('Observation Point','interpreter','latex')
ylabel('$ |\int_0^T \frac {\partial G}{\partial m} dt| $','interpreter','latex', 'fontsize',20)
ylim([1e-10 1e-5])

str3 = sprintf('Results/SOAB1D5/Plots/D5_case_%d_sigma_%d_dG_dm_integrated.eps', case_num, k);
print(str3, '-depsc2')
str5 = sprintf('Results/SOAB1D5/Plots/D5_case_%d_sigma_%d_dG_dm_integrated.fig', case_num, k);
savefig(str5)

    end
end

%% Figure 6: Plots for dG/dm integrated over time
clear variables
for case_num = 1
    
if case_num ==1
    test = [0.01, 0.05, 0.07, 0.1, 0.15, 0.17, 0.2, 0.25,0.27, 0.3];
elseif case_num == 2
    test = [0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2];
elseif case_num == 3
    test = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1];
end

    for k = [1 5 10]
str1 =sprintf('Results/SOAB1D7/DA Results/SOAB1D7_case_%d_amplitude_%d.mat', case_num, k);
str2 =sprintf('Results/SOAB1D7/SOA Results/SOAB1D7_results_case_%d_amplitude_%d.mat', case_num, k);
load(str1)
load(str2)

N = length(beta_exct0);
cfl = 1/3;
xmax = 3; 
xmin = -xmax;
tmax = 6;
tmin = 0;
delta_x = (xmax - xmin)/N;
dt = delta_x*cfl;
T = tmin:dt:tmax;


sens_at_obs = eta_hat(x0_inds,:);
sens_int_T = trapz(T, sens_at_obs, 2);

figure(k); 
semilogy(abs(sens_int_T), 'linewidth',2)
set(gca,'FontSize',14);
ylim([1e-11 1e-6])


xlabel('Observation Point','interpreter','latex')
ylabel('$| \int_0^T \frac {\partial G}{\partial m} dt| $','interpreter','latex', 'fontsize',20)


str3 = sprintf('Results/SOAB1D7/Plots/case_%d_amplitude_%d_dG_dm_integrated.eps', case_num, k);
print(str3, '-depsc2')
str5 = sprintf('Results/SOAB1D7/Plots/case_%d_amplitude_%d_dG_dm_integrated.fig', case_num, k);
savefig(str5)

    end
end

%% Figure 7: Plots for dG/dm integrated over time
clear variables
for case_num = 2
    
if case_num ==1
    test = [0.01, 0.05, 0.07, 0.1, 0.15, 0.17, 0.2, 0.25,0.27, 0.3];
elseif case_num == 2
    test = [0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2];
elseif case_num == 3
    test = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1];
end

    for k = [1 4 9]
str1 =sprintf('Results/SOAB1D7/DA Results/SOAB1D7_case_%d_amplitude_%d.mat', case_num, k);
str2 =sprintf('Results/SOAB1D7/SOA Results/SOAB1D7_results_case_%d_amplitude_%d.mat', case_num, k);
load(str1)
load(str2)

N = length(beta_exct0);
cfl = 1/3;
xmax = 3; 
xmin = -xmax;
tmax = 6;
tmin = 0;
delta_x = (xmax - xmin)/N;
dt = delta_x*cfl;
T = tmin:dt:tmax;


sens_at_obs = eta_hat(x0_inds,:);
sens_int_T = trapz(T, sens_at_obs, 2);

figure(k); 
semilogy(abs(sens_int_T), 'linewidth',2)
set(gca,'FontSize',14);
ylim([1e-9 1e-2])


xlabel('Observation Point','interpreter','latex')
ylabel('$| \int_0^T \frac {\partial G}{\partial m} dt| $','interpreter','latex', 'fontsize',20)


str3 = sprintf('Results/SOAB1D7/Plots/case_%d_amplitude_%d_dG_dm_integrated.eps', case_num, k);
print(str3, '-depsc2')
str5 = sprintf('Results/SOAB1D7/Plots/case_%d_amplitude_%d_dG_dm_integrated.fig', case_num, k);
savefig(str5)

    end
end

