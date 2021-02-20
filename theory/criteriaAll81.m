% This script finds all the relevant criteria in Extended Data Table 1:
% Effect of cost is independent of partner density and the two populations 
% have separate carrying capacities. 

% This script performs the following tasks:
% 1) solve the steady state solutions of the 27 models.
% 2) 27 figures are generated with each figure number corresponding to 
%    its model number
%      a) for the models that violate asusmptions, demonstrate how one or
%         more assumpations are violated
%      b) for models that satisfy the assumptions, solve for the criteria
%         and validate with vector fields and time course simulations

% The four assumptions of model behavior: 
% a) Benefit shall be positively dependent on partner density.
% b) Cost shall decrease growth rate or carrying capacity.
% c) Stress shall produce negative growth of populations with some 
%    parameter combinations. 
% d) Negative growth of a population shall be potentially counteracted by
%    benefit produced by its partner but further strengthened by cost. 

% The corresponding symbols for parameters in Table 1 are as follows: 
% d is the symbol for delta (stress)
% b is the symbol for beta (benefit)
% c is the symbol for epsilon (cost)
% b2 is the symbol for beta prime (a second parameter for benefit)

% updated 2/19/18
% Copyright 2018, Feilun Wu, all rights reserved. 

%% set up the symbols
syms x1 x2
syms d b c b2 z

% Eliminate the cases when parameters equal to 0 to simplify the analyses
assume(d>0); assume(b>0); assume(c>0); assume(b2>0);

addpath('.\I_CriteriaTables\baseForm_27')
addpath('.\I_CriteriaTables\shareCarryingCap_27')
addpath('.\I_CriteriaTables\densityDependentCost_27')


%% Solve for all steady state solutions
ssAll = cell(81,1);
for i = 1:81
    ssAll{i} = eval(['solve(ODE' sprintf('%02d',i) '(0,[x1;x2],[d,b,c,b2])==[0;0],[x1,x2]);']);
end

% Verify that the model solutions are symmetric
% Needs to check symmetry of model if the steady state solutions are not
% symmetric. Note that some symmetric models can yield asymmetric
% solutions, such as ODE28. This can be due to the existence of line steady
% states. 
sanityCheck = zeros(81,1);
for i = 1:81
    sanityCheck(i) = isequal(sort(ssAll{i}.x1),sort(ssAll{i}.x2));
end

%% Showcases of how a model violates assumptions
% The figures of this section contains three panels
% 1) The top panel is the dx1/dt vs x2. This plot should be monotonically 
% increasing to satisfy the model assmuption.
% 2) The middle panel is the dx1/dt vs c (cost). This plot should be monotonically 
% decreasing to satisfy the model assmuption.
% 3) The bottom panel is an instance of growth curves. Because of symmetry,
% the growth curves of the two populations overlap. This growth curve
% should be bounded. 

% The plots generated in this section demonstrate that these models violate 
% at least one of above requirements. 

modelViol_base = [2 3 5 6 8 9 10 12 13 14 16]; % model structures that violate assumptions
modelViol = [modelViol_base modelViol_base+27 modelViol_base+54];

cTurnover1 = [8, 16, 35, 43, modelViol_base+54];

for i = 1:length(modelViol)
    if ismember(modelViol(i),cTurnover1)
        % These model formulations have cost implemented as linear death
        % rate and is greater than 0. 
        par = [1.5 2 0.1 1];
    else
        
        par = [1.5 2 1.1 1];
    end

    eval(['assmptExp(@ODE' sprintf('%02d',modelViol(i)) ',par);']);
end

%% For other models, find the expressions of criteria through investigating steady state solutions

modelSatf_ = [1 4 7 11 15 17:27]; % model structures that satisfy assumptions
modelSatf = [modelSatf_ modelSatf_+27 modelSatf_+54];
Bfunction = cell(81,1);
ssHigh = cell(81,1);
ind_needFurtherInv = [0 0]; tmpCount = 0;
cTurnover2 = [7, 17, 34, 44, modelSatf_+54];
undSqrtExp = cell(81,1);

for i = 1:length(modelSatf)
    MI = modelSatf(i);
    
    if ismember(modelSatf(i),cTurnover2)
        % These model formulations have cost implemented as linear death
        % rate and is greater than 0. 
        par = [1.2 200 0.1 2];
    else
        
        par = [1.2 200 1.1 2];
    end
    
    ssValues = getSS(ssAll{MI},[d,b,c,b2],par);

    n = length(ssValues);
    if n == 0
        ind_needFurtherInv(tmpCount+1) = MI;
        tmpCount = tmpCount+1;
    else
        isEql=(ssValues(:,1)==ssValues(:,2));
        [~,I]=sort(real(ssValues(:,1).*isEql),'descend');
        ssHigh{MI} = ssAll{MI}.x1(I(1));
        undSqrtExp{MI} = underSqrt(ssHigh{MI});
        if ~isempty(undSqrtExp{MI}{1})&& contains(undSqrtExp{MI}{1},'d')
            eval(['Bfunction{MI} = solve((',undSqrtExp{MI}{1},')==0,d);']);
            if length(Bfunction{MI})==2  
                if double(subs(Bfunction{MI}(1) - Bfunction{MI}(2),[d,b,c,b2],[1.2,2,0.1,1]))<0
                    Bfunction{MI} = (Bfunction{MI}(1));
                else
                    Bfunction{MI} = (Bfunction{MI}(2));
                end
            end
        elseif strfind(char(ssHigh{MI}),'root(')==1
            % a function that calculates the discriminant of a polynomial
            discrim = findsDiscriminant(ssHigh{MI});
            Bfunction{MI} = solve(discrim==0,d);
        else 
            ind_needFurtherInv(tmpCount+1) = MI;
            tmpCount = tmpCount+1;
        end
    end
end


%% Check the models that could not be solved using expressions under the square root
% Thes two model's non-trivial stable steady states are constant. However, the unstable
% steady state that determines the separatrix of the system is determined by parameters.
% For these two systems to coexist, the initial density has to be greater than the
% unstable steady state. 
syms x0
%%====================ODE01====================
par =[1.2 5 1.1 1]; % a test parameter set
d_t = par(1); b_t = par(2); c_t = par(3);  % unwrap "par"
Bfunction{1} = (x0*b+1)/(c);
x0_crit = (c_t*d_t - 1)/b_t; % calculate the boundary x0

figure(1)
subplot(3,1,1)
plotVF_SS(@ODE01,ssAll{1},par,[0 2 0 2])

subplot(3,1,2)
x0_v = x0_crit*1.1;
runModel(@ODE01,[x0_v,x0_v],par,200);
title('meets criterion: B/\delta>1')
xlabel('Time (\tau)');ylabel('Density')

subplot(3,1,3)
x0_v = x0_crit*0.9;
runModel(@ODE01,[x0_v,x0_v],par,200);
title('violates criterion: B/\delta<1')
xlabel('Time (\tau)');ylabel('Density')

%%====================ODE04====================
par =[1.2 5 1.1 1]; % a test parameter set
d_t = par(1); b_t = par(2); c_t = par(3); % unwrap "par"
Bfunction{4} = (x0_v*b+1);
x0_crit = (d_t - 1)/b_t; % calculate the boundary x0

figure(4)
subplot(3,1,1)
plotVF_SS(@ODE04,ssAll{4},par,[0 2 0 2])

subplot(3,1,2)
x0_v = x0_crit*1.1;
runModel(@ODE04,[x0_v,x0_v],par,200);
title('meets criterion: B/\delta>1')
xlabel('Time (\tau)');ylabel('Density')

subplot(3,1,3)
x0_v = x0_crit*0.9;
runModel(@ODE04,[x0_v,x0_v],par,200);
title('violates criterion: B/\delta<1')
xlabel('Time (\tau)');ylabel('Density')

%%====================ODE28====================
par =[1.2 5 1.1 1]; % a test parameter set
d_t = par(1); b_t = par(2); c_t = par(3);  % unwrap "par"
Bfunction{28} = (x0_v*b+1)/(c);
x0_crit = (c_t*d_t - 1)/b_t; % calculate the boundary x0

figure(28)
subplot(3,1,1)
plotVF_SS(@ODE28,ssAll{28},par,[0 2 0 2])

subplot(3,1,2)
x0_v = x0_crit*1.1;
runModel(@ODE28,[x0_v,x0_v],par,200);
title('meets criterion: B/\delta>1')
xlabel('Time (\tau)');ylabel('Density')

subplot(3,1,3)
x0_v = x0_crit*0.9;
runModel(@ODE28,[x0_v,x0_v],par,200);
title('violates criterion: B/\delta<1')
xlabel('Time (\tau)');ylabel('Density')

%%====================ODE31====================
par =[1.2 5 1.1 1]; % a test parameter set
d_t = par(1); b_t = par(2); c_t = par(3); % unwrap "par"
Bfunction{31} = (x0_v*b+1);
x0_crit = (d_t - 1)/b_t; % calculate the boundary x0

figure(31)
subplot(3,1,1)
plotVF_SS(@ODE31,ssAll{31},par,[0 2 0 2])

subplot(3,1,2)
x0_v = x0_crit*1.1;
runModel(@ODE31,[x0_v,x0_v],par,200);
title('meets criterion: B/\delta>1')
xlabel('Time (\tau)');ylabel('Density')

subplot(3,1,3)
x0_v = x0_crit*0.9;
runModel(@ODE31,[x0_v,x0_v],par,200);
title('violates criterion: B/\delta<1')
xlabel('Time (\tau)');ylabel('Density')

%%====================ODE55====================
par =[1.2 5 1.1 1]; % a test parameter set
d_t = par(1); b_t = par(2); c_t = par(3);  % unwrap "par"
Bfunction{55} = (x0_v*b+1)/(c*x0_v+1);
x0_crit = (d_t - 1)/(b_t-d_t*c_t); % calculate the boundary x0

figure(55)
subplot(3,1,1)
plotVF_SS(@ODE55,ssAll{55},par,[0 2 0 2])

subplot(3,1,2)
x0_v = x0_crit*1.1;
runModel(@ODE55,[x0_v,x0_v],par,200);
title('violates criterion: B/\delta>1')
xlabel('Time (\tau)');ylabel('Density')

subplot(3,1,3)
x0_v = x0_crit*0.9;
runModel(@ODE55,[x0_v,x0_v],par,200);
title('meets criterion: B/\delta<1')
xlabel('Time (\tau)');ylabel('Density')

%%====================ODE58====================
par =[1.2 5 1.1 1]; % a test parameter set
d_t = par(1); b_t = par(2); c_t = par(3); % unwrap "par"
Bfunction{58} = (x0_v*b+1);
x0_crit = (d_t - 1)/b_t; % calculate the boundary x0

figure(58)
subplot(3,1,1)
plotVF_SS(@ODE58,ssAll{58},par,[0 2 0 2])

subplot(3,1,2)
x0_v = x0_crit*1.1;
runModel(@ODE58,[x0_v,x0_v],par,200);
title('meets criterion: B/\delta>1')
xlabel('Time (\tau)');ylabel('Density')

subplot(3,1,3)
x0_v = x0_crit*0.9;
runModel(@ODE58,[x0_v,x0_v],par,200);
title('violates criterion: B/\delta<1')
xlabel('Time (\tau)');ylabel('Density')

%% Test the model criteria

countS = 0;
rootOfExp = [0 0];

for i = 1:length(modelSatf)
    MI = modelSatf(i); % model index
    if ismember(MI,cTurnover2)
        % These model formulations have cost implemented as linear death
        % rate and is greater than 0. 
        par = [1.5 3 0.1 1.2];
    else  
        par = [1.5 3 1.1 1.2];
    end
    
    if ~ismember(MI,ind_needFurtherInv)
        if length(Bfunction{MI})==1
            dcrit = double(subs(Bfunction{MI},[d,b,c,b2],par));
            eval(['testCriterion(@ODE' sprintf('%02d',MI) ',ssAll{MI},dcrit,par)']);
        elseif length(Bfunction{MI})>=1
            dcrits = double(subs(Bfunction{MI},[d,b,c,b2],par));
            posReal = abs(imag(dcrits))<eps; % find the real roots         
            [sortedDel,I1] = sort(dcrits.*posReal,'ascend');
            smallestPos = find(sortedDel>0,1);
            delCrit = dcrits(I1(smallestPos));

            eval(['testCriterion(@ODE' sprintf('%02d',MI) ',ssAll{MI},delCrit,par)']);
            Bfunction{MI} = Bfunction{MI}(I1(smallestPos));
            countS = countS+1;
            rootOfExp(countS) = MI;
        end
    end
end

%% check whether all 81 models are investigated
size( get(0,'Children'))

%% plot how Bfunction change with b and c
% All top right corners have higher benefit

critInd = ~cellfun(@isempty,Bfunction);
totCrit = sum(sum(critInd,2))-length(rootOfExp);
n = 20;
[b_mesh, c1_mesh] = meshgrid(linspace(2,4,n),linspace(1,2,n));
[~, c2_mesh] = meshgrid(linspace(2,4,n),linspace(0.1,0.2,n));
ind = 1;
figure(100)
for i = 1:81
    if ismember(MI,cTurnover2)
        % These model formulations have cost implemented as linear death
        % rate and is greater than 0. 
        c_mesh = c2_mesh;
    else  
        c_mesh = c1_mesh;
    end
    if critInd(i)==1 && ~ismember(i,rootOfExp) && ~ismember(i,ind_needFurtherInv)
        fH = matlabFunction(Bfunction{i},'Vars',[d b c b2]);
        subplot(6,7,ind)
        imagesc(fH(1.2,b_mesh,c_mesh,1));
        set(gca,'xTick',[1 n],'xTickLabel',[2 4],'yTick',[1 n],'yTickLabel',[min(c_mesh(:)) max(c_mesh(:))])
        title(i)
        xlabel('\beta')
        ylabel('\epsilon')
        ind = ind+1;
    end
end

