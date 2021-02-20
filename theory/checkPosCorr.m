syms d b c b2

tspan = 300;
dcrit = double(subs(Bfunction{MI},[d,b,c,b2],par));

delta_crit=ssAll{72};

% testCriterion(model,modelss,delta_crit,par,modelInd)

figure()

subplot(2,2,1)
d_this = 1.1*1.05;

% par = [08 3 0.1 1.2];
%%
par(4) = 1.8;
[~,y] = runModel(@ODE72,[0.5,0.5],par,tspan);
title('violates criterion: B/\delta<1')
xlabel('Time (\tau)');ylabel('Density')
axis([0 tspan 0 1.5])

%%
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