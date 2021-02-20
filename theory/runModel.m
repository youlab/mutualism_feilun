function [t,y]=runModel(ODE,y0,params,tend)
% This function runs a time-course simulation of an ODE model and plot the
% time course. 
% 
% INPUTS
% ODE: the dydt function of the ODE model.
% y0: initial conditions
% params: model parameters
% tend: end time point of simulation


options=odeset('MaxStep',1,'NonNegative',[1,2]);
[t,y]=ode45(ODE,[0 tend],y0,options,params);
plot(t,y)

end