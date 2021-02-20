function VectorF(func,params,axis_custom)
% This function plots 2D vector field of a ODE model.

% INPUTS
% func: the dydt function
% params: model parameters
% axis_custom: the axis range of the vector field

x1 = linspace(0.00001,axis_custom(2),30);
x2 = linspace(0.00001,axis_custom(4),30);
[x(1,:,:),x(2,:,:)] = meshgrid(x1,x2);

dx=func(0,x,params);
cosa=squeeze(dx(1,:,:));
sina=squeeze(dx(2,:,:));
length=sqrt(cosa.^2+sina.^2);

quiver(x1,x2,cosa./length,sina./length);
xlabel('x1')
ylabel('x2')
axis('square');
axis(axis_custom)



%% one time course
% options=odeset('MaxStep',1,'NonNegative',[1,2]);
% [sx,sy] = ode15s(@I_ODE,[0,100],[0.1;0.2],options,[0.5,100,0.5],1);
% plot(sx,sy)
% plot(sy(:,1),sy(:,2),'r');
% plot(0.001,0.001,'ko')
% 
% axis('square');
% axis([0 1 0 1])