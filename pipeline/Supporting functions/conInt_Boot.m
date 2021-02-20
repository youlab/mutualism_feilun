function [Biass,Variances,Biass2,Variances2, RSD, STD, Amp] = conInt_Boot(X,Y,n,grid,ker,kerV,BoxC,B_real)

% This function performs the bootstrapping for a single set of kernel 
% parameter and caculates how different the calibrated B is from the real 
% B. It is called by getIntervals() to loop through the bootstrapping 
% process for a set of kernel parameters.

% X, Y are the input data
% n is the number of bootstrap to perform
% grid is the v (system variable) mesh
% ker is the kernel function name
% kerV is the kernel only includes v but not delta
% BoxC is box constraint used in SVM training
% B_real is the true B

B = zeros([length(grid),n]);

for i = 1:n
    inds = datasample(1:length(Y),ceil(length(Y)*0.9),'Replace',true);
    modelTemp = fitcsvm(X(inds,:),Y(inds,:),'KernelFunction',ker,'BoxConstraint',BoxC);
    B(:,i) = getB(X(inds,:),Y(inds,:),modelTemp,grid,kerV,X,Y); % adding X and Y is only for deciding the directionality
end

Biass = mean((B-repmat(B_real(:),1,n)).^2,2); % bias
B_m = mean(B,2); % mean surface of all bootstrapped B
Variances= mean((B-repmat(B_m(:),1,n)).^2,2); % variance of all bootstrapped B
Amp = max(B_m)-min(B_m); % amplitude of the mean B surface

STD = std(B,[],2); % standard deviation of the bootstrapped B
RSD = STD./abs(B_m(:)); % relative standard deviation of the bootstrapped B

B_normm = zeros([length(grid),n]); % normalize the n-column B matrix column by column
for i = 1:n
    r = polyfit(B(:,i),B_real(:),1);
    B_normm(:,i) = (B(:,i)*r(1))+r(2);
end

% Bias, mean and variance after rescaling the bootstrapped B to the real B
Biass2 = mean((B_normm-repmat(B_real(:),1,n)).^2,2);
B_m2 = mean(B_normm,2);
Variances2= mean((B_normm-repmat(B_m2(:),1,n)).^2,2);

end