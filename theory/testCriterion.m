function testCriterion(model,modelss,delta_crit,par,modelInd)

% This function first finds the boundary delta value based on Bfunction and
% then plot two cases. To verify the correctness of the criterion, when
% B/d<1, both the time-series simulation and vector field should only
% indicate a single steady state at (0,0). When B/d>1, the two populations
% should both be positive and the vector field contains a stable steady 
% state where both populations have positive densities. 

% INPUTS
% model: the ODE model of interest
% modelss: the steady state solution of the ODE model
% Bfunction: the theoretical B function of b, c, and b2. It represents
% effective benefit.
% par: a specific set of model parameter values. par is in the order of 
% [d, b, c, b2]

% Outputs:
% This function plots two scenarios: when stress (d) is greater (the top 
% row) and less (the bottom row) than B. It plots the time course (left
% panels) and corresponding vector field with steady states indicated by
% black circles. 

syms d b c b2

tspan = 300;

if nargin==4
    modelInd = regexp(func2str(model),'\d*','Match');
    figure(str2double(modelInd{1}))
elseif nargin==5
    figure(modelInd)
end

subplot(2,2,1)
d_this = delta_crit*1.05;
par(1) = d_this;
[~,y] = runModel(model,[0.5,0.5],par,tspan);
title('violates criterion: B/\delta<1')
xlabel('Time (\tau)');ylabel('Density')
axis([0 tspan 0 1.5])

subplot(2,2,2)
plotVF_SS(model,modelss,par,[0 1 0 1])
plot(y(:,1),y(:,2))

subplot(2,2,3)
d_this = delta_crit*0.95;
par(1) = d_this;
[~,y] = runModel(model,[1,1],par,tspan); % use high initial density to eliminate the effect of saddle points
title('meets criterion: B/\delta>1')
xlabel('Time (\tau)');ylabel('Density')
axis([0 tspan 0 1.5])

subplot(2,2,4)
plotVF_SS(model,modelss,par,[0 1 0 1])
plot(y(:,1),y(:,2))