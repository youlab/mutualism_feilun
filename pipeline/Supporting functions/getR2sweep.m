function [rsq_,rsq_2,rsq_3,a1,a2]=getR2sweep(k1_,k2_,X_norm,Y,ker,BoxC,B_real,sigmaX,miuX)   

% This function go through all the kernel parameter combinations and
% calculate the R2 values of the trained B surface compared withe the real
% B. 

% k1_, k2_ are kernel parameters of interest
% X_norm, Y are input data
% ker is the kernel function name
% BoxC is box constraint used in SVM training
% B_real is the true B
% sigmaX, miuX are mean and standard deviation of X


global k1 k2

rsq_ = zeros(length(k1_),length(k2_));
rsq_2 = zeros(length(k1_),length(k2_));
rsq_3 = zeros(length(k1_),length(k2_));
a1 = zeros(length(k1_),length(k2_));
a2 = zeros(length(k1_),length(k2_));

for i = 1:length(k1_)
    for j = 1:length(k2_)
        k1 = k1_(i);
        k2 = k2_(j);
        SVMmodel = fitcsvm(X_norm,Y,'KernelFunction',ker,'BoxConstraint',BoxC);
        eval(['B = getB(X_norm,Y,SVMmodel,X_norm(:,1:2),@' ker 'V,X_norm,Y);']);
     
        B_scaledback = ((B.*sigmaX)+miuX);
        delta_back = ((X_norm(:,3).*sigmaX)+miuX);

        rsq_(i,j) = 1 - sum((B_real(:) - B_scaledback(:)).^2)/sum((B_real(:) - mean(B_real(:))).^2); % directly compare B and Breal
        
        r = polyfit(B_scaledback,B_real(:),1);
        B_linFit = B_scaledback.*r(1)+r(2);
        rsq_2(i,j)= 1 - sum((B_real(:) - B_linFit).^2)/sum((B_real(:) - mean(B_real(:))).^2); % compare a fitted B with Breal (proportionality)
        a1(i,j) = r(1);
        
        XXX=B_scaledback(:)./delta_back(:)-1;
        YYY=B_real(:)./delta_back(:)-1;
        rBd = polyfit(XXX(:),YYY(:),1);
        rsq_3(i,j)= 1 - sum((YYY(:) - XXX.*rBd(1)-rBd(2)).^2)/sum((YYY(:) - mean(YYY(:))).^2); % compare Breal/d and B/d
        a2(i,j) = rBd(1);

    end
end