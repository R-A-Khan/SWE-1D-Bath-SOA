
load('resampling.mat')
load('YB_sobol.mat')


% Find indices of all undefined values
ind = find(isnan(Y_beta));
for i = 1: length(ind)
disp(['k = ', sprintf('%i', ind(i)),', Bath amplitude = ', sprintf('%0.3f', amp_beta(ind(i))),', IC amplitude = ', sprintf('%0.4f', amp_ic(ind(i))),', ph_shift = ', sprintf('%0.4f', ph_shift(ind(i))),', Y_beta = ', sprintf('%0.4e', Y_beta(ind(i))),', Y_eta = ', sprintf('%0.4e', Y_eta(ind(i)))])
end

% Randomly generate indices to remove from Y
% match up sample size to smallest N base rate such that no undefined vals
% in YA, YB, YC

for k = 1:7
    ind2(k) = randi(2993);
end

% Remove indices for undefined values and randomly generated indices
YB_beta_new = Y_beta;
YB_beta_new(ind) = [];
YB_beta_new(ind2) = [];

YB_eta_new = Y_eta;
YB_eta_new(ind) = [];
YB_eta_new(ind2) = [];

XB_new = XB;
XB_new(ind,:) = [];
XB_new(ind2,:) = [];


% New values should all contain no undefined wavles, YA and YB should be 
% N x 1, YC should be M*N X 1, corresponding to parameters in
% XA (N x 1), XB (N x 1), XC(M*N x 1)


save('YB_sobol_new.mat','YB_beta_new', 'YB_eta_new', 'XB_new');