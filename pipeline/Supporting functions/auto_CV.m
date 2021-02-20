function [classLoss] = auto_CV(X,Y,k1_,k2_,kerName,BoxC)

% This script finds the mean CV loss and std of CV loss for all parameter
% combinations

% X, Y are input data
% k1_ and k2_ are kernel parameters of interest
% kerName is the kernel function
% BoxC is box constraint used in SVM training

global k1 k2

classLoss = zeros(10,length(k1_),length(k2_)); % 10 is the default CV fold number

for i = 1:length(k1_)
    for j = 1:length(k2_)
        k1 = k1_(i);
        k2  = k2_(j);
        SVMModel = fitcsvm(X, Y, 'KernelFunction', kerName,'BoxConstraint',BoxC);
        CVSVMModel = crossval(SVMModel);
        classLoss(:,i,j) = kfoldLoss(CVSVMModel,'mode','individual');
    end
end

