function [meanLoss_top,params,inds]=getTop(meanLoss,Var,lamb,k1,k2,n)
% this script finds the parameters for training the top-scored calibrated B

% meanLoss: the mean CV loss
% Var: the variance across all bootstrapped B
% lamb: the parameter that dictates the values of the scoring metric
% k1, k2: the kernel parameters
% n: the number of top parameters of interest

[meanLoss_s, ind]= sort((meanLoss(:)*lamb)+ (1-lamb)*Var(:));
meanLoss_top = meanLoss_s(1:n);

k1_m = repmat(k1',[1,length(k2)]);
k2_m = repmat(k2,[length(k1),1]);
k1_l = k1_m(:);
k2_l = k2_m(:);

% if meanLoss and Var are values from many kernels, repeat the linearized k
% (k1_l and k2_l) again by the number of kernels (nn1/nn2).
nn1 = length(meanLoss);
nn2 = length(k1_l);
k1_l = repmat(k1_l(:),nn1/nn2,1);
k2_l = repmat(k2_l(:),nn1/nn2,1);

inds=ind(1:n);
params = [k1_l(inds) k2_l(inds)];

end