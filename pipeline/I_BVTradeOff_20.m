% This script trains the 20 models shown in Extended Data Figure 4b. 

% updated 2/19/18
% Copyright 2018, Feilun Wu, all rights reserved. 

%% setup the basic parameters
numKP = 20;
numModel = 20;
k1_ = logspace(-2,2,numKP);
k2_ = logspace(-2,2,numKP);
BoxC = 10;
ker = {'myLin2','myQuad2','myCub2','mySigm2'};
global k1 k2

colors =  [        
    0         0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

addpath('./1 Extended data fig 4')
addpath('./Supporting functions')
addpath('./fourKernels')

load establishRationale20.mat
%% An example that shows how to train the B surface using a particular kernel 
% and a particular set of model parameters. 
% Here, no kernel parameter screening is performed. 

k1 =1;
k2 =1;

miuX = mean(X{5});
sigmaX = std(X{5});
X_norm=(X{5}-miuX)./sigmaX;

modelTemp = fitcsvm(X_norm,Y{5},'KernelFunction',ker{1},'BoxConstraint',BoxC);
Bthis= getB(X_norm,Y{5},modelTemp,X_norm(:,1:2),@myLin2V);

figure(1)
subplot(1,4,1)
imagesc(reshape(Bthis,10,10));colorbar
title('B_{calibrated}')
subplot(1,4,2)
imagesc(reshape(B_real{5},10,10));colorbar
title('B_{real}')
subplot(1,4,3)
imagesc(reshape(X_norm(:,3)*sigmaX(3)+miuX(3),10,10));colorbar
title('\delta')
subplot(1,4,4)
imagesc(reshape(Y{5},10,10));colorbar
title('Y')

%% training the SVM models for all 20 data sets

% ***if want to directly use the pre-ran results, do not run this section, 
% instead: load update_getIntervals.mat and 
% update_getIntervals_cvloss_R2.mat directly ***

load update_getIntervals.mat
load update_getIntervals_cvloss_R2.mat

% ** The following portion may run more than 30 hours on a personal
% computer **

% % initialize all the metrics of a B surface
% R2_ = cell(numModel,4);
% R2_loose = cell(numModel,4);
% R2_Bd = cell(numModel,4);
% classLoss_ = cell(numModel,4);
% Variances_= cell(numModel,4);
% Biass_= cell(numModel,4);
% Variances2_= cell(numModel,4);
% Biass2_= cell(numModel,4);
% RSD_= cell(numModel,4);
% STD_= cell(numModel,4);
% Amp_= cell(numModel,4);
% a1_= cell(numModel,4);
% a2_= cell(numModel,4);
% 
% for i = 1:20
%     X_this = X{i};
%     Y_this = Y{i};
%     Breal_this = B_real{i};
% 
%     miuX = mean(X_this);
%     sigmaX = std(X_this);
%     X_norm=(X_this-miuX)./sigmaX;
% 
%         for j = 1:4
%             tic
%             ker_this = ker{j};
%             [classLoss_{i,j}] = auto_CV(X_norm,Y_this,k1_,k2_,ker_this,1);
%             [Biass_{i,j},Variances_{i,j},Biass2_{i,j},Variances2_{i,j}, RSD_{i,j},STD_{i,j},Amp_{i,j}] = getIntervals(X_norm,Y_this,k1_,k2_,ker_this,X_norm(:,1:2),BoxC,500,Breal_this);
%             [R2_{i,j},R2_loose{i,j},R2_Bd{i,j},a1_{i,j},a2_{i,j}]= getR2sweep(k1_,k2_,X_norm,Y_this,ker_this,BoxC,Breal_this,sigmaX(3),miuX(3)); 
%             toc
%         end 
%     disp(i)
% end
% 
% % save the calibration results
% save update_getIntervals.mat Biass_ Variances_ Biass2_ Variances2_ RSD_ STD_ Amp_
% save update_getIntervals_cvloss_R2.mat classLoss_ R2_ R2_loose R2_Bd a1_ a2_

%% get single values for each set of kernel parameter

varT = zeros(numKP,numKP);
biasT = zeros(numKP,numKP);
var2T = zeros(numKP,numKP);
bias2T = zeros(numKP,numKP);
STDT = zeros(numKP,numKP);
RSDT = zeros(numKP,numKP);

var_ = cell(numModel, 4);
bias_ = cell(numModel, 4);
var2_ = cell(numModel, 4);
bias2_ = cell(numModel, 4);
rsd_ = cell(numModel, 4);
std_ = cell(numModel, 4);
cvloss_ = cell(numModel, 4);
cvlossStd_ = cell(numModel, 4);

for h = 1:20
    for l = 1:4
        for i = 1:numKP
            for j = 1:numKP
                varT(i,j)= mean(Variances_{h,l}{i,j});
                biasT(i,j)= mean((Biass_{h,l}{i,j}));
                var2T(i,j)= mean(Variances2_{h,l}{i,j});
                bias2T(i,j)= mean((Biass2_{h,l}{i,j}));
                RSDT(i,j)= mean(RSD_{h,l}{i,j});
                STDT(i,j)= mean((STD_{h,l}{i,j}));
            end
        end
        var_{h,l}=varT;
        bias_{h,l}=biasT;
        var2_{h,l}=var2T;
        bias2_{h,l}=bias2T;
        rsd_{h,l}=RSDT;
        std_{h,l}=STDT; 
        cvloss_{h,l} = squeeze(mean(classLoss_{h,l}));
        cvlossStd_{h,l} = squeeze(std(classLoss_{h,l}));
    end
end

%% Visualize all the 20 sets of data
% The titles indicate the number of the models.
% A model contains 3 plots, from left to right: true B landscape; delta;
% classification of coexistence or collapse

figure(2)
ind = 1;

for i = 1:20
    subplot(4,15,ind)
    imagesc(B_real{i})
%     colorbar
    title(i)
    subplot(4,15,ind+1)
    imagesc(reshape(X{i}(:,3),10,10))
%     colorbar
    subplot(4,15,ind+2)
    imagesc(reshape(Y{i},10,10))
    ind = ind+3;
end

%% Plot how average R2 changes with lambda
% this plot helps identify the optimal lambda to use for ranking the B
% surfaces. 

nn = 1000;
ind = 1;
lambda = linspace(0,1,nn);
mmm = 3; % get the top mmm calibrated B
R2save = zeros(nn,mmm,80);
OptR2 = zeros(nn,20,mmm);

% numKP = 20;

for i = 1:20
    cvloss4 = zeros(numKP*numKP*4,1);
    var4 = zeros(numKP*numKP*4,1);
    R24 = zeros(numKP*numKP*4,1);
    
    for j = 1:4 % this step ensures that simpler kernels are prioritized in the ranking process
        cvloss4((1+(j-1)*numKP^2):(j*numKP^2),1) = cvloss_{i,j}(:);
        var4((1+(j-1)*numKP^2):(j*numKP^2),1) = var_{i,j}(:);
        R24((1+(j-1)*numKP^2):(j*numKP^2),1) = R2_loose{i,j}(:);
    end
    
    for k = 1:nn
        [~,~,I] = getTop(cvloss4,var4,lambda(k),k1_,k2_,mmm);
        OptR2(k,i,:)= (R24(I)); 
    end

    ind = ind+1;
end

tempOptR2 = zeros(nn,mmm*20);
for k = 1:nn
    tt= OptR2(k,:,:);
    tempOptR2(k,:) = tt(:);
end

figure(3)
hold on
plot(lambda,nanmean(tempOptR2,2))
set(gca,'yscale','log')
axis([0 1 0.65 0.73])
xlabel('\lambda')
ylabel('R^2')

%% plot supplemental figure (all four kernels together)
lamOptAll = 0.8; % the optimal lambda obtained from above. 
Ind = 1;
XXX=cell(4,4);
YYY=cell(4,4);
percR2 = zeros(20,4);
minScore = zeros(20,4); % find min score for all kernels and models

ind = 1;
for i = 1:20
    for j = 1:4
        figure(4)
        subplot(5,4,i)
        hold on
        rr = R2_loose{i,j};
        varTemps = var_{i,j};

        XXX{i,j}= lamOptAll*cvloss_{i,j}(:)+(1-lamOptAll)*varTemps(:);
        YYY{i,j}= rr(:);
        minScore(i,j) = min(XXX{i,j});
        [~,I]=sort(XXX{i,j},'ascend');
        percR2(i,j) = nanmean(R2_loose{i,j}(I(1:1)));
        scatter(XXX{i,j},YYY{i,j},'o','filled','MarkerFaceAlpha',3/8,'MarkerFaceColor',colors(j,:))
        axis([1E-5 100 0.0 1])
        set(gca,'xscale','log')
        ind = ind+1;

    end
        hold off
end

%% find all top B for each kernel with lamda=0.8 
% (80 calibraed B compared with 20 real B)
% Column 1 and 6 are real B surfaces.
% The real B surface is followed by 4 calibrated B with
% different kernels. From left to right, the four kernels are linear,
% quadratic, cubic and sigmoidal. 

lamOpt = 0.8;
topR2 = zeros(5,4);
figure(5)

ind = 1;
for i = 1:20
    X_this = X{i};
    Y_this = Y{i};
    Breal_this = B_real{i};

    subplot(10,10,ind)
    imagesc(Breal_this)
    
    miuX = mean(X_this);
    sigmaX = std(X_this);
    X_norm=(X_this-miuX)./sigmaX;
    ind = ind +1;
    
    for j = 1:4
        ker_this = ker{j};
        
        [~,kthis,~] = getTop(cvloss_{i,j}(:),var_{i,j}(:),lamOpt,k1_,k2_,1);
        k1 = kthis(1);
        k2 = kthis(2);
        
        modelTemp = fitcsvm(X_norm,Y_this,'KernelFunction',ker_this,'BoxConstraint',BoxC);
        eval(['Btemp = getB(X_norm,Y_this,modelTemp,X_norm(:,1:2),@' ker_this 'V,X_norm,Y_this);']);        
        
        B_cali = ((Btemp.*sigmaX(3))+miuX(3)); 
        
        r = polyfit(B_cali,Breal_this(:),1); 
        B_linFit = B_cali.*r(1)+r(2);
        topR2(i,j)= 1 - sum((B_real{i}(:) - B_linFit).^2)/sum((B_real{i}(:) - mean(B_real{i}(:))).^2); % compare a fitted B with Breal (proportionality)

        subplot(10,10,ind)
        imagesc(reshape(Btemp,10,10))
        
        ind = ind +1;
    end
end
%% find all top 4 B with all kernels considered with lamda=0.8 
% (80 calibraed B compared with 20 real B)
% Column 1 and 6 are real B surfaces.
% The real B surface is followed by 4 calibrated B with
% different kernels. From left to right, the four kernels are linear,
% quadratic, cubic and sigmoidal. 

% ** Note that although a majority of calibrated B capture the shape of true
% B, a few calibrated B do not share similar profile with the real
% B, which indicates that the calibration process can benefit from further
% improvements. **

lamOpt = 0.8;
topR2 = zeros(5,4);
figure(6)
mmm = 4;

ind = 1;
for i = 1:20
    X_this = X{i};
    Y_this = Y{i};
    Breal_this = B_real{i};
    
    subplot(10,10,ind)
    imagesc(Breal_this)
    
    miuX = mean(X_this);
    sigmaX = std(X_this);
    X_norm=(X_this-miuX)./sigmaX;
    ind = ind +1;
    
    for j = 1:4 % this step ensures that simpler kernels are prioritized in the ranking process
        cvloss4((1+(j-1)*numKP^2):(j*numKP^2),1) = cvloss_{i,j}(:);
        var4((1+(j-1)*numKP^2):(j*numKP^2),1) = var_{i,j}(:);
        R24((1+(j-1)*numKP^2):(j*numKP^2),1) = R2_loose{i,j}(:);
    end
    [~,kthis,I] = getTop(cvloss4,var4,lamOpt,k1_,k2_,mmm);

    for j = 1:mmm
        ker_this = ker{ceil(I(j)/(numKP.^2))};
        k1 = kthis(j,1);
        k2 = kthis(j,2);
        
        modelTemp = fitcsvm(X_norm,Y_this,'KernelFunction',ker_this,'BoxConstraint',BoxC);
        eval(['Btemp = getB(X_norm,Y_this,modelTemp,X_norm(:,1:2),@' ker_this 'V,X_norm,Y_this);']);        
        
        B_cali = ((Btemp.*sigmaX(3))+miuX(3));
        
        r = polyfit(B_cali,Breal_this(:),1);
        B_linFit = B_cali.*r(1)+r(2);
        topR2(i,j)= 1 - sum((B_real{i}(:) - B_linFit).^2)/sum((B_real{i}(:) - mean(B_real{i}(:))).^2); % compare a fitted B with Breal (proportionality)

        subplot(10,10,ind)
        imagesc(reshape(B_cali,10,10))
        
        ind = ind +1;
    end
end

