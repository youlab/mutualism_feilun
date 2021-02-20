% This script finds the calibrated B for the QS-based mutualism system

% updated 2/19/18
% Copyright 2018, Feilun Wu, all rights reserved. 

%% Calibration setup
% set up the kernel parameters
numKP = 20;
numModel = 1;
k1_ = logspace(-2,2,numKP);
k2_ = logspace(-2,2,numKP);

global k1 k2

BoxC = 10;
ker = {'myLin2','myQuad2','myCub2','mySigm2'};
addpath('./fourKernels')
addpath('./Supporting functions')
addpath('./2 M1M2')
load M1M2_96.mat

% normalize the input data
miuX = mean(X);
sigmaX = std(X);
X_norm = (X-miuX)./sigmaX;

colors =  [        
    0         0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

%% Train SVM with 4 different kernels

load M1M2_results.mat
% classLoss_ = cell(1,4);
% Variances_ = cell(1,4);
% 
% % ** the following code may run more than 4 hours on a personal computer
% for j = 1:4
%     tic
%     ker_this = ker{j};
%     [classLoss_{j}] = auto_CV(X_norm,Y,k1_,k2_,ker_this,1);
%     [Variances_{j},~,~,~,~] = getIntervals_use(X_norm,Y,k1_,k2_,ker_this,X_norm(:,1:2),BoxC,500);
%     toc
% end
% save M1M2_results.mat classLoss_ Variances_

%% calculate a single variance for a set of parameter

varT = zeros(numKP,numKP); % temporary variable

% initialize variables
var_ = cell(numModel, 4);
cvloss_ = cell(numModel, 4);
cvlossStd_ = cell(numModel, 4);

for h = 1:numModel
    for l = 1:4
        for i = 1:numKP
            for j = 1:numKP
                varT(i,j)= mean(Variances_{h,l}{i,j});
            end
        end

        var_{h,l}=varT;
        cvloss_{h,l} = squeeze(mean(classLoss_{h,l}));
        cvlossStd_{h,l} = squeeze(std(classLoss_{h,l}));
    end
end

%% plot the top 5 

lamOpt = 0.8; 
mmm = 5; % the number of top B of interest

% set up the labels 
labelsIPTG = reshape(X(:,1),8,12);labelsATC = reshape(X(:,2),8,12);
ticksPosIPTG = 1:5:12;ticksPosATC = 1:3:8;

% integrate the results from 4 kernels
% this step ensures that simpler kernels are prioritized in the ranking process
for j = 1:4 
    cvloss4((1+(j-1)*numKP^2):(j*numKP^2),1) = cvloss_{j}(:);
    var4((1+(j-1)*numKP^2):(j*numKP^2),1) = var_{j}(:);
end
[~,kthis,I] = getTop(cvloss4,var4,lamOpt,k1_,k2_,mmm);

% plot the 5 calibrated B and their RSD
figure(101)
for j = 1:mmm
    ker_this = ker{ceil(I(j)/(numKP.^2))};
    k1 = kthis(j,1);
    k2 = kthis(j,2);

    modelTemp = fitcsvm(X_norm,Y,'KernelFunction',ker_this,'BoxConstraint',BoxC);
    eval(['Btemp = getB(X_norm,Y,modelTemp,X_norm(:,1:2),@' ker_this 'V,X_norm,Y);']);        

    B_cali = ((Btemp.*sigmaX(3))+miuX(3));

    subplot(2,mmm,j)
    imagesc(reshape(B_cali,8,12))
    set(gca,'xTick',ticksPosIPTG,'xTickLabel',labelsIPTG(1,ticksPosIPTG))
    set(gca,'yTick',ticksPosATC,'yTickLabel',labelsATC(ticksPosATC,1))
    colorbar
    
    eval(['[B_u,B_l,B_m,Variances,RSD,~,~] = conInt_Boot_use(X_norm,Y,10,[X_norm(:,[1 2])],ker_this,@',ker_this,'V,BoxC);'])

    subplot(2,mmm,j+mmm)

    imagesc(reshape(RSD,8,12))
    set(gca,'xTick',ticksPosIPTG,'xTickLabel',labelsIPTG(1,ticksPosIPTG))
    set(gca,'yTick',ticksPosATC,'yTickLabel',labelsATC(ticksPosATC,1))
    caxis([0 1])
    colorbar
    
    figure(123)
    subplot(1,mmm,j)
    plot(B_cali(:)./X(:,3),yy(:),'o')
    axis([0.8 1.3 0 0.5])
end

% add the titles of the two rows
figure(101)
subplot(2,mmm,1)
ylabel('B_{cali}')
subplot(2,mmm,1+mmm)
ylabel('RSD')

%%
