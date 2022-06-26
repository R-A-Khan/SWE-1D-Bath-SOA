%% Plots for HV = F

addpath(genpath('/Users/Ramsha/Dropbox/McMaster University/PhD Research/Matlab Code/SOA Implementation/Results'))

%%
clear all
case_num = 1;
 test =[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];

for k = 1:length(test)
str1 =sprintf('Results/SOAB1D5/DA Results/SOAB1D5_case_%d_sigma_%d.mat', case_num, k);
str2 =sprintf('Results/SOAB1D5/SOA Results/SOAB1D5_results_case_%d_sigma_%d.mat', case_num, k);
load(str1)
load(str2)

figure(k); 
plot(Hv,'r'); hold on; plot(F, 'b');hold off;
str3 = sprintf('Bathymetry Gaussian with Standard Deviation %0.2f',test(k));
title(str3)
xlabel('x','interpreter','latex')
legend('Hv','F','location','northeast')
str4 = sprintf('Results/SOAB1D5/Plots/Case_%d_Hv_F_sigma_%d.eps',case_num, k);
print(str4, '-depsc2')
str5 = sprintf('Results/SOAB1D5/Plots/Case_%d_Hv_F_sigma_%d.fig',case_num, k);
savefig(str5)
end

%% Plots for dG/dM


% addpath(genpath('/Users/Ramsha/Dropbox/McMaster University/PhD Research/Matlab Code/SOA Implementation/Results'))
clear all
case_num = 1;
 test =[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
for k = 1:length(test)
str1 =sprintf('Results/SOAB1D5/DA Results/SOAB1D5_case_%d_sigma_%d.mat', case_num, k);
str2 =sprintf('Results/SOAB1D5/SOA Results/SOAB1D5_results_case_%d_sigma_%d.mat', case_num, k);
load(str1)
load(str2)
figure(k); 

subplot(3,1,1)
str4 = sprintf('Sensitivity: Bathymetry Gaussian with Standard Deviation %0.2f',test(k));
plot(eta_hat(x0_inds(1),:),'r'); 
legend('Obs 1','location','northeast')
xlabel('time step','interpreter','latex')
title(str4, 'interpreter','latex')

subplot(3,1,2)
plot(eta_hat(x0_inds(round(n_obs/2)),:),'b');
str_plot1 = sprintf('0bs %d',round(n_obs/2)); 
legend( str_plot1,'location','northeast')
xlabel('time step','interpreter','latex')

subplot(3,1,3)
plot(eta_hat(x0_inds(end),:),'k'); 
str_plot2 = sprintf('0bs %d',n_obs); 
legend(str_plot2,'location','northeast')
xlabel('time step','interpreter','latex')

str3 = sprintf('Results/SOAB1D5/Plots/Case_%d_dG_dm_sigma_%d.eps',case_num, k);
print(str3, '-depsc2')
str5 = sprintf('Results/SOAB1D5/Plots/Case_%d_dG_dm_sigma_%d.fig',case_num, k);
savefig(str5)
end