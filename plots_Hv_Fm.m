%% Plots for HV = F

addpath(genpath('/Users/Ramsha/Dropbox/McMaster University/PhD Research/Matlab Code/SOA Implementation/Results'))

%%
clear all
for case_num = 1:3;

str1 =sprintf('Results/SOAB1D1/DA Results/Case_DAB1D21_case_%d.mat', case_num);
str2 =sprintf('Results/SOAB1D1/SOA Results/results_bicgstabl_case_%d.mat', case_num);
load(str1)
load(str2)

figure(case_num); 
plot(Hv,'r'); hold on; plot(F, 'b');hold off;
% str3 = sprintf('Bathymetry Gaussian with Standard Deviation %0.2f',test(k));
% title(str3)
xlabel('x','interpreter','latex')
legend('Hv','F','location','northeast')
str4 = sprintf('Results/SOAB1D1/Plots/Case_%d_Hv_F.eps',case_num);
print(str4, '-depsc2')
str5 = sprintf('Results/SOAB1D1/Plots/Case_%d_Hv_F.fig',case_num);
savefig(str5)
end


%% Plots for dG/dM


% addpath(genpath('/Users/Ramsha/Dropbox/McMaster University/PhD Research/Matlab Code/SOA Implementation/Results'))
clear all
for case_num = 1:3;

str1 =sprintf('Results/SOAB1D1/DA Results/Case_DAB1D21_case_%d.mat', case_num);
str2 =sprintf('Results/SOAB1D1/SOA Results/results_bicgstabl_case_%d.mat', case_num);
load(str1)
load(str2)
figure(case_num); 

subplot(3,1,1)
% str4 = sprintf('Sensitivity: Bathymetry Gaussian with Standard Deviation %0.2f',test(k));
plot(eta_hat(x0_inds(1),:),'r'); 
legend('Obs $1$','location','northeast', 'interpreter','latex')
xlabel('time step','interpreter','latex')
% title(str4, 'interpreter','latex')

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

str3 = sprintf('Results/SOAB1D1/Plots/Case_%d_dG_dm.eps',case_num);
print(str3, '-depsc2')
str5 = sprintf('Results/SOAB1D1/Plots/Case_%d_dG_dm.fig',case_num);
savefig(str5)
end

%%

for case_num = 1:3

str1 =sprintf('Results/SOAB1D1/DA Results/Case_DAB1D21_case_%d.mat', case_num);
str2 =sprintf('Results/SOAB1D1/SOA Results/results_bicgstabl_case_%d.mat', case_num);
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
figure(case_num); 
plot(sens_int_T, 'linewidth',2)
set(gca,'FontSize',14);

xlabel('Observation Point','interpreter','latex')
ylabel('$ \int_0^T \frac {\partial G}{\partial m} dt $','interpreter','latex', 'fontsize',20)


str3 = sprintf('Results/SOAB1D1/Plots/Case_%d_dG_dm_integrated.eps',case_num);
print(str3, '-depsc2')
str5 = sprintf('Results/SOAB1D1/Plots/Case_%d_dG_dm_integrated.fig',case_num);
savefig(str5)
end
