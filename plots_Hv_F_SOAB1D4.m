%% Plots for HV = F

addpath(genpath('/Users/Ramsha/Dropbox/McMaster University/PhD Research/Matlab Code/SOA Implementation/Results'))

%%
clear all
case_num = 1;
 test =[ 5 10 15 25 35 45 75];

for k = 1:length(test)
str1 =sprintf('Results/SOAB1D4/DA Results/SOAB1D4_case_%d_nObs_%d.mat', case_num,test(k));
str2 =sprintf('Results/SOAB1D4/SOA Results/SOAB1D4_results_case_%d_obs_no_%d.mat', case_num, k);
load(str1)
load(str2)

figure(k); 
plot(Hv,'r'); hold on; plot(F, 'b');hold off;
str3 = sprintf(' %d Observation Points',n_obs);
title(str3)
xlabel('x','interpreter','latex')
legend('Hv','F','location','northeast')

% str4 = sprintf('Results/SOAB1D4/Plots/Case_%d_Hv_F_case_%d.eps',case_num, test(k));
% print(str4, '-depsc2')

str5 = sprintf('Results/SOAB1D4/Plots/Case_%d_Hv_F_case_%d.fig',case_num, test(k));
savefig(str5)
end

%% Plots for dG/dM


% addpath(genpath('/Users/Ramsha/Dropbox/McMaster University/PhD Research/Matlab Code/SOA Implementation/Results'))
clear all
case_num = 1;
test =[ 5 10 15 25 35 45 75];
for k = 1:length(test)
str1 =sprintf('Results/SOAB1D4/DA Results/SOAB1D4_case_%d_nObs_%d.mat', case_num,test(k));
str2 =sprintf('Results/SOAB1D4/SOA Results/SOAB1D4_results_case_%d_obs_no_%d.mat', case_num, k);
load(str1)
load(str2)
figure(k); 

subplot(3,1,1)
str4 = sprintf('Sensitivity with %d Observation Points',n_obs); 
plot(eta_hat(x0_inds(1),:),'r'); 
legend('Obs 1','location','northeast')
xlabel('time step','interpreter','latex')
title(str4, 'interpreter','latex')

subplot(3,1,2)
plot(eta_hat(x0_inds(round(n_obs/2)),:),'b');
legend( '0bs 23','location','northeast')
xlabel('time step','interpreter','latex')

subplot(3,1,3)
plot(eta_hat(x0_inds(end),:),'k'); 
legend('Obs 45','location','northeast')
xlabel('time step','interpreter','latex')

% str3 = sprintf('Results/SOAB1D4/Plots/Case_%d_dG_dm_case_%d.eps',case_num, test(k));
% print(str3, '-depsc2')

str5 = sprintf('Results/SOAB1D4/Plots/Case_%d_dG_dm_case_%d.fig',case_num, test(k));
savefig(str5)
end


%% Plots for dG/dm integrated over time
clear all
for case_num = 1:3

    
test =[ 5 10 15 25 35 45 75];


    for k = 1:length(test)
str1 =sprintf('Results/SOAB1D4/DA Results/SOAB1D4_case_%d_nObs_%d.mat', case_num,test(k));
str2 =sprintf('Results/SOAB1D4/SOA Results/SOAB1D4_results_case_%d_obs_no_%d.mat', case_num, k);
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

xlabel('Observation Point','interpreter','latex')
ylabel('$ \int_0^T \frac {\partial G}{\partial m} dt $','interpreter','latex', 'fontsize',20)

str4 = sprintf('%d Observation Points',n_obs);
title(str4, 'interpreter', 'latex')

str3 = sprintf('Results/SOAB1D4/Plots/case_%d_N_obs_%d_dG_dm_integrated.eps', case_num, test(k));
print(str3, '-depsc2')
str5 = sprintf('Results/SOAB1D4/Plots/case_%d_N_obs_%d_dG_dm_integrated.fig', case_num, test(k));
savefig(str5)

    end
end

%% Plots for obs
clear all
for case_num = 1:3
test =[ 5 10 15 25 35 45 75];

for k = 1:length(test)
str1 =sprintf('Results/SOAB1D4/DA Results/SOAB1D4_case_%d_nObs_%d.mat', case_num,test(k));
load(str1)
figure(k);
plot(X,eta_exct0(:,1), 'k-','LineWidth',1.5); hold on
plot(X, 0.01*beta_exct0(:,1,1)-0.01, 'k--','LineWidth',1.5);
scatter(X(x0_inds), Y_ex(x0_inds), 'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7],'LineWidth',1.5); hold off
% ylim([-2e-2 4e-3])
xlim([min(X) max(X)])
set(gca,'YTick', [])
% title('Bathymetry and IC (not to scale)', 'fontsize',8, 'fontweight', 'bold' )
set(gca,'FontSize',18)
str4 = sprintf('%d Observation Points',n_obs);
title(str4, 'interpreter', 'latex', 'fontsize', 14)

str3 = sprintf('Results/SOAB1D4/Plots/case_%d_N_obs_%d_obs_number.eps', case_num, test(k));
print(str3, '-depsc2')
str5 = sprintf('Results/SOAB1D4/Plots/case_%d_N_obs_%d_obs_number.fig', case_num, test(k));
savefig(str5)
end
end
