%% Plots for HV = F

%%
clear all
case_num = 1;
 test =[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];

for k = 1:length(test)
str1 =sprintf('Results/SOAB1D6/DA Results/SOAB1D6_case_%d_sigma_%d.mat', case_num, k);
str2 =sprintf('Results/SOAB1D6/SOA Results/SOAB1D6_results_case_%d_sigma_%d.mat', case_num, k);
load(str1)
load(str2)


figure(k); 
plot(Hv,'r'); hold on; plot(F, 'b');hold off;
str3 = sprintf('Bathymetry Gaussian with Standard Deviation %0.2f',test(k));
title(str3)
xlabel('x','interpreter','latex')
legend('Hv','F','location','northeast')
str4 = sprintf('Results/SOAB1D6/Plots/Case_%d_Hv_F_sigma_%d.eps',case_num, k);
print(str4, '-depsc2')
str5 = sprintf('Results/SOAB1D6/Plots/Case_%d_Hv_F_sigma_%d.fig',case_num, k);
savefig(str5)
end

%% Plots for dG/dM


% addpath(genpath('/Users/Ramsha/Dropbox/McMaster University/PhD Research/Matlab Code/SOA Implementation/Results'))
clear all
case_num = 1;
 test =[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
for k = 1:length(test)
str1 =sprintf('Results/SOAB1D6/DA Results/SOAB1D6_case_%d_sigma_%d.mat', case_num, k);
str2 =sprintf('Results/SOAB1D6/SOA Results/SOAB1D6_results_case_%d_sigma_%d.mat', case_num, k);
load(str1)
load(str2)
figure(k); 

subplot(3,1,1)
str4 = sprintf('Sensitivity with Shifted IC: \n Bathymetry Gaussian with Standard Deviation %0.2f',test(k));
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

str3 = sprintf('Results/SOAB1D6/Plots/Case_%d_dG_dm_sigma_%d.eps',case_num, k);
print(str3, '-depsc2')
str5 = sprintf('Results/SOAB1D6/Plots/Case_%d_dG_dm_sigma_%d.fig',case_num, k);
savefig(str5)
end


%% Plots for dG/dm integrated over time
clear variables
for case_num = 1

    
 test =[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
    for k = 1:length(test)
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
plot(sens_int_T, 'linewidth',2)
set(gca,'FontSize',14);
str4 = sprintf('Bathymetry Gaussian with Standard Deviation %0.2f',test(k));
title(str4, 'interpreter', 'latex')


xlabel('Observation Point','interpreter','latex')
ylabel('$ \int_0^T \frac {\partial G}{\partial m} dt $','interpreter','latex', 'fontsize',20)


str3 = sprintf('Results/SOAB1D6/Plots/case_%d_sigma_%d_dG_dm_integrated.eps', case_num, k);
print(str3, '-depsc2')
str5 = sprintf('Results/SOAB1D6/Plots/case_%d_sigma_%d_dG_dm_integrated.fig', case_num, k);
savefig(str5)

    end
end

%% Plots for obs
clear variables
for case_num = 1
 test =[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
    for k = 1:length(test)
str1 =sprintf('Results/SOAB1D6/DA Results/SOAB1D6_case_%d_sigma_%d.mat', case_num, k);
load(str1)

figure(k);
plot(X,eta_exct0(:,1), 'k-','LineWidth',1.5); hold on
plot(X, 0.01*beta_exct0(:,1,1)-0.01, 'k--','LineWidth',1.5);
% scatter(X(x0_inds), Y_ex(x0_inds), 'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7],'LineWidth',1.5); hold off
% ylim([-2e-2 4e-3])
xlim([min(X) max(X)])
set(gca,'YTick', [])

set(gca,'FontSize',18)
str4 = sprintf('Bathymetry Standard Deviation %0.2f',test(k));
title(str4, 'interpreter', 'latex', 'fontsize', 14)

str3 = sprintf('Results/SOAB1D6/Plots/case_%d_N_obs_%d_bath_shape.eps', case_num, k);
print(str3, '-depsc2')
str5 = sprintf('Results/SOAB1D6/Plots/case_%d_N_obs_%d_bath_shape.fig', case_num, k);
savefig(str5)
    end
end
