%% Fixing initial plots for HV = F
for k = 1:3
test = subplot(3,1,k);
test.XTick = [0 384 769];
test.XTickLabel = {'0','$\frac{T}{2}$','$T$'}
test.TickLabelInterpreter = 'latex';
subplot(3,1,k);xlabel('$t$', 'interpreter', 'latex') 
set(gca,'FontSize',16);

end

savefig('Case_3_dG_dm.fig')
print('Case_3_dG_dm.eps', '-depsc2')



%% Fixing plots for Case I
open('Case_3_dG_dm_integrated.fig')
set(gca,'FontSize',20);
axis square
savefig('Case_3_dG_dm_integrated.fig')
print('Case_3_dG_dm_integrated.eps', '-depsc2')

%% Fixing plots for Question 2 Case I

case_num = 3 ;

set(0,'DefaultFigureWindowStyle','docked')
pts = [1 10 20];
for k = pts
str = sprintf('case_%d_x0_%d_dG_dm_integrated.fig', case_num, k);
open(str)
set(gca,'FontSize',20);
title('')
axis square
savefig(str)
str2 = sprintf('case_%d_x0_%d_dG_dm_integrated.eps', case_num, k);
print(str2, '-depsc2')
end

%%
set(0,'DefaultFigureWindowStyle','docked')
case_num = 3 ;
pts = [1 10 20];
for k = pts
str = sprintf('case_%d_x0_%d_obs_position.fig', case_num, k);
open(str)
set(gca,'FontSize',20);
xlabel('$x$', 'interpreter', 'latex') 
axis square
savefig(str)
str2 = sprintf('case_%d_x0_%d_obs_position.eps', case_num, k);
print(str2, '-depsc2')
end


%% Case 2c (SOAB1D6)

case_num = 1 ;

set(0,'DefaultFigureWindowStyle','docked')
pts = [2 5 8];
for k = pts
str = sprintf('case_%d_sigma_%d_dG_dm_integrated.fig', case_num, k);
open(str)
set(gca,'FontSize',20);
title('')
axis square
savefig(str)
str2 = sprintf('case_%d_sigma_%d_dG_dm_integrated.eps', case_num, k);
print(str2, '-depsc2')
end

%%
set(0,'DefaultFigureWindowStyle','docked')
case_num = 1 ;
pts = [2 5 8];
for k = pts
str = sprintf('case_%d_N_obs_%d_bath_shape.fig', case_num, k);
open(str)
set(gca,'FontSize',20);
xlabel('$x$', 'interpreter', 'latex') 
axis square
title('')
savefig(str)
str2 = sprintf('case_%d_sigma_%d_bath_shape.eps', case_num, k);
print(str2, '-depsc2')
end

%% New plots for bath shape in SOAB1D6
set(0,'DefaultFigureWindowStyle','docked')

for case_num = 1:3
test = [2 5 8];
    for k = test  
str1 =sprintf('Results/SOAB1D6/DA Results/SOAB1D6_case_%d_sigma_%d.mat', case_num, k);
load(str1)
figure;
plot(X,eta_exct0(:,1), 'k-','LineWidth',1.5); hold on
plot(X, 0.01*beta_exct0(:,1,1)-0.01, 'k--','LineWidth',1.5);
% scatter(X(x0_inds), Y_ex(x0_inds), 'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7],'LineWidth',1.5); hold off
% ylim([-2e-2 4e-3])
xlim([min(X) max(X)])
set(gca,'YTick', [])
xlabel('$x$', 'interpreter', 'latex') 
% title('Bathymetry and IC (not to scale)', 'fontsize',8, 'fontweight', 'bold' )
set(gca,'FontSize',20)
axis square

str3 = sprintf('Results/SOAB1D6/Plots/case_%d_sigma_%d_bath_shape.eps', case_num, k);
print(str3, '-depsc2')
str5 = sprintf('Results/SOAB1D6/Plots/case_%d_sigma_%d_bath_shape.fig', case_num, k);
savefig(str5)
    end
end

%% Plots for 2ci (SOAB1D5)

case_num = 1 ;

set(0,'DefaultFigureWindowStyle','docked')
pts = [2 5 8];
for k = pts
str = sprintf('case_%d_N_obs_%d_bath_shape.fig', case_num, k);
open(str)
set(gca,'FontSize',20);
xlabel('$x$', 'interpreter', 'latex') 
axis square
title('')
savefig(str)
str2 = sprintf('D5_case_%d_N_obs_%d_bath_shape.eps', case_num, k);
print(str2, '-depsc2')
end

%% SOAB1D7

%% Case 2c (SOAB1D7)

case_num = 3 ;

set(0,'DefaultFigureWindowStyle','docked')
pts = [1 5 10];
for k = pts
str = sprintf('case_%d_amplitude_%d_dG_dm_integrated.fig', case_num, k);
open(str)
set(gca,'FontSize',20);
title('')
axis square
savefig(str)
str2 = sprintf('case_%d_amplitude_%d_dG_dm_integrated', case_num, k);
print(str2, '-depsc2')
end

%%
set(0,'DefaultFigureWindowStyle','docked')
case_num = 2 ;
pts = [1 4 9];
for k = pts
str = sprintf('case_%d_amplitude_%d_bath_shape.fig', case_num, k);
open(str)
set(gca,'FontSize',20);
xlabel('$x$', 'interpreter', 'latex') 
axis square
title('')
savefig(str)
str2 = sprintf('case_%d_amplitude_%d_bath_shape', case_num, k);
print(str2, '-depsc2')
end



