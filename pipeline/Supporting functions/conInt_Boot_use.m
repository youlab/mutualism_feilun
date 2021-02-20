function [B_u,B_l,B_m,Variances,RSD,B,P]=conInt_Boot_use(X,Y,n,grid,ker,kerV,BoxC)

% This function is a simple version of conInt_Boot() and do not make
% comparison between the calibrated B and real B. 

% It performs the bootstrapping for a single set of kernel 
% parameter. It is called by getIntervals_use.m to loop through the 
% bootstrapping process for a set of kernel parameters.

% X, Y are the input data
% n is the number of bootstrap to perform
% grid is the v (system variable) mesh
% ker is the kernel function name
% kerV is the kernel only includes v but not delta
% BoxC is box constraint used in SVM training

B = zeros([length(grid),n]);

for i = 1:n
    inds = datasample(1:length(Y),ceil(length(Y)),'Replace',true);
    modelTemp = fitcsvm(X(inds,:),Y(inds,:),'KernelFunction',ker,'BoxConstraint',BoxC);
    B(:,i) = getB(X(inds,:),Y(inds,:),modelTemp,grid,kerV,X,Y); % adding X and Y is only for deciding the directionality
end

B_m = mean(B,2);

Variances = mean((B-repmat(B_m(:),1,n)).^2,2);

STD = std(B,[],2);
RSD = STD./B_m;
SEM = STD/sqrt(n);               % Standard Error

ts = tinv([0.025  0.975],n-1);      % T-Score

B_u = mean(B,2) + ts(2).*SEM;                      % Confidence Intervals
B_l = mean(B,2) + ts(1).*SEM; 

P = abs(B_u-B_l)./B_m; % 95% confidence interval compared with mean

end