%% Results with smoothing for DA



load('Results/SOAB1D1/DA Results/Case_DAB1D21_case_1.mat');
err_filt_1 = err;
cost_filt_1 = cost;
beta_opt_filt_1 = beta_optimum;


load('Results/SOAB1D1/DA Results/Case_DAB1D21_case_2.mat');
err_filt_2 = err;
cost_filt_2 = cost;
beta_opt_filt_2 = beta_optimum;



load('Results/SOAB1D1/DA Results/Case_DAB1D21_case_3.mat');
err_filt_3 = err;
cost_filt_3 = cost;
beta_opt_filt_3 = beta_optimum;


iter = 1:iter_max;
skip=iter_max/20;
marker1 = {'k-+', 'k-o', 'k-*', 'k-x', 'k-s', 'k-d', 'k-^', 'k-v', 'k->', 'k-<', '-kp', 'k-h'};
marker2 = {'k--+','k--o','k--*','k--x','k--s','k--d','k--^','k--v','k-->','k--<','k--p','k--h'};

figure(1);
% semilogy(iter(1:skip:end), err_no_filt_iv(1:skip:end,1), marker1{1},'markersize',10,'linewidth',1.5); grid on; hold on;
semilogy(iter(1:skip:end), err_filt_3(1:skip:end,1), marker1{1},'markersize',10,'color','b','linewidth',1.5);
set(gca,'FontSize',18); grid on;
xlabel('Iteration n', 'interpreter', 'latex'); ylabel('Relative $L^2$ error for $\beta(x)$','interpreter','latex');
xlim([0 iter_max]);
% ylim([1e-2 1e0]);
axis square
savefig('Results/SOAB1D1/Plots/err_case3.fig')
print -depsc2 Results/SOAB1D1/Plots/err_case3.eps
