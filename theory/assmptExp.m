function assmptExp(model,par)

% This function tests whether the model satisfy assumption d) and whether
% it generates unbounded growth.

modelInd = regexp(func2str(model),'\d*','Match');

figure(str2double(modelInd{1}))
subplot(3,1,1) % part 1 of assumption d)
x2_ = linspace(0,1,1000);
x1_ = ones(1,1000)/2;
dydt_ = model(0,[x1_;x2_],par);
plot(x2_,dydt_(1,:))
xlabel('x_2')
ylabel('dx_1/d\tau')
title(func2str(model))

subplot(3,1,2) % part 2 of assumption d)
c_ = linspace(1,2,1000);
dydt2_ = zeros(1000,2);
for i = 1:1000
    par(3) = c_(i);
    dydt2_(i,:) = model(0,[0.2;0.2],par);
end
plot(c_,dydt2_(:,1))
xlabel('\epsilon')
ylabel('dx_1/d\tau')

subplot(3,1,3) % test unbounded growth
runModel(model,[1,1],par,100);
xlabel('Time (\tau)')
ylabel('Density')