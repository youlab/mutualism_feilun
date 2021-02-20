function plotVF_SS(func,ss,par,axesRange)

% This function plots the vector field together with steady states. 
% The steady states, regardless of stability, are plotted as open circles.

% INPUTS
% func: the ODE model of interest
% ss: symbolic expression of steady state solutions of the ODE model
% par: a specific set of model parameter values. par is in the order of 
% [d, b, c, b2]
% axesRange: the axis range of the vector field 

syms x1 x2
syms d b c b2

ss_val = getSS(ss,[d b c b2],par); % each row is one steady state

hold on
VectorF(func,par,axesRange)

for i = 1:length(ss_val)
    if sum(ss_val(i,:)>=0)==2 && sum(abs(imag(ss_val(i,:)))<eps)==2
        plot(real(ss_val(i,1)),real(ss_val(i,2)),'ko')
    end
end
