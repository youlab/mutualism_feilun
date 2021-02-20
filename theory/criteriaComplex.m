% Calculates and verify the criteria with increasing complexity

% This script 
% a) finds all the relevant criteria in Extended Data Table 2:
%    models with relaxed assumptions
% b) verifies the criteria by using time-course simualtions and vector
%    field

% updated 2/19/18
% Copyright 2018, Feilun Wu, all rights reserved. 

%% set up the symbols and simulation options
syms x1 x2 x0 x n
syms d1 d2 b1 b2 c1 c2 a1 a2 d01 d02 p
syms d b c a d0

% Eliminate the cases when parameters equal to 0 to simplify the analyses
assume(d>0); assume(b>0); assume(c>0); assume(a>0); assume(x>0);
assume(in(n,'integer') & n>2)

addpath('.\II_CriteriaComplexity')

tspan = 300;

%% Add competition
% a: parameter that dictates the strength of competition. 

ss_comp_82 = solve(ODE_base_comp(0,[x1;x2],[d d b b c c a a 0 0 1])==[0;0],[x1,x2]);
par82s1 = [1.2 200 1.1 0.8];
par82s2 = [1.2 5 1.1 0.5];
par82_long = [1.2 1.2 5 5 1.1 1.1 0.5 0.5 0 0 1];
parSym82 = [d b c a];

% find the high steady state
ssValues82 = getSS(ss_comp_82,parSym82,par82s1);
isEql=(ssValues82(:,1)==ssValues82(:,2));
[~,I]=sort(real(ssValues82(:,1).*isEql),'descend');
ssHigh82 = ss_comp_82.x1(I(1));
% get the effective benefit function through investigating the expression
% under the square root 
undSqrtExp82 = underSqrt(ssHigh82);
%=============B function================
eval(['Bfunction82 = solve((',undSqrtExp82{1},')==0,d);']);
%=======================================

% verify the criterion
dcrit = double(subs(Bfunction82,parSym82,par82s2));
figure(82)
subplot(2,2,1)
d_this = dcrit*1.05;
par82_long(1) = d_this; par82_long(2) = d_this;
[~,y] = runModel(@ODE_base_comp,[0.5,0.5],par82_long,tspan);
title('violates criterion: B/\delta<1')
xlabel('Time (\tau)');ylabel('Density')
axis([0 tspan 0 1])

subplot(2,2,2)
hold on
VectorF(@ODE_base_comp,par82_long,[0 1 0 1])
par82s2(1) = d_this;
ss_val = getSS(ss_comp_82,parSym82,par82s2); 
for i = 1:length(ss_val)
    if sum(ss_val(i,:)>=0)==2 && sum(abs(imag(ss_val(i,:)))<eps)==2
        plot(real(ss_val(i,1)),real(ss_val(i,2)),'ko')
    end
end
plot(y(:,1),y(:,2))

subplot(2,2,3)
d_this = dcrit*0.95;
par82_long(1) = d_this; par82_long(2) = d_this;
[~,y] = runModel(@ODE_base_comp,[0.5,0.5],par82_long,tspan);
title('meets criterion: B/\delta>1')
xlabel('Time (\tau)');ylabel('Density')
axis([0 tspan 0 1])

subplot(2,2,4)
hold on
VectorF(@ODE_base_comp,par82_long,[0 1 0 1])
par82s2(1) = d_this;
ss_val = getSS(ss_comp_82,parSym82,par82s2); 
for i = 1:length(ss_val)
    if sum(ss_val(i,:)>=0)==2 && sum(abs(imag(ss_val(i,:)))<eps)==2
        plot(real(ss_val(i,1)),real(ss_val(i,2)),'ko')
    end
end
plot(y(:,1),y(:,2))

%% A criterion that considers initial density

ss_comp_83 = solve(ODE_base_comp(0,[x1;x2],[d d b b c c 0 0 0 0 1])==[0;0],[x1,x2]);
par83s1 = [1.2 200 1.1];
par83s2 = [1.2 5 1.1];
par83_long = [1.2 1.2 5 5 1.1 1.1 0 0 0 0 1];
parSym83 = [d b c];

% identifies the steady state that dictates the saddle point
ssValues83 = getSS(ss_comp_83,parSym83,par83s1);
isEql=(ssValues83(:,1)==ssValues83(:,2));
[~,I]=sort(real(ssValues83(:,1).*isEql),'descend');
ssSaddle83 = ss_comp_83.x1(I(2)); 
%=============B function================
Bfunction83 = solve(ssSaddle83==x0,d);
%=======================================

% verifies the criterion
x0crit =  double(subs(ssSaddle83,[d,b,c],par83s2));
figure(83)
subplot(3,1,2)
[~,y1] = runModel(@ODE_base_comp,[x0crit*1.1,x0crit*1.1],par83_long,tspan);
subplot(3,1,3)
[~,y2] = runModel(@ODE_base_comp,[x0crit*0.9,x0crit*0.9],par83_long,tspan);
subplot(3,1,1)
hold on
VectorF(@ODE_base_comp,par83_long,[0 1 0 1])

ss_val = getSS(ss_comp_83,parSym83,par83s2); 
for i = 1:length(ss_val)
    if sum(ss_val(i,:)>=0)==2 && sum(abs(imag(ss_val(i,:)))<eps)==2
        plot(real(ss_val(i,1)),real(ss_val(i,2)),'ko')
    end
end
plot(y1(:,1),y1(:,2)) % the trace that starts below the saddle point
plot(y2(:,1),y2(:,2)) % the trace that starts above the saddle point

%% Include asymmetry 
% stress (d1, d2), benefit (b1, b2), cost (c1, c2) and ratio of growth rate
% (p) can all represent system asymmetry. 

parSym84 = [d1 d2 b1 b2 c1 c2 1 1 0 0 p];
ss_comp_84 = solve(ODE_base_comp(0,[x1;x2],parSym84)==[0;0],[x1,x2]);
par84_long =  [1.2 1.2 100 100 1.1 1.1 1 1 0 0 1];
par84_long2 = [1.2 1.2 5 5 1.1 1.1 1 1 0 0 1.2];

% find the accurate criterion
ssValues84 = getSS(ss_comp_84,parSym84,par84_long);
isEql=(ssValues84(:,1)==ssValues84(:,2));
[~,I]=sort(real(ssValues84(:,1).*isEql),'descend');
ssHigh84 = ss_comp_84.x1(I(1));
undSqrtExp84 = underSqrt(ssHigh84);
%=============B function================
eval(['Bfunction84 = solve(' , undSqrtExp84{1}, '==0,d1);']);
%=======================================

% an approximation of the criterion
theta1 = b2/b1; theta2 = (c2*d2)/(c1*d1); % parameters indicating asymmetry
myFunc = p*(b1+1/theta1+1)^2/(4*b1*c1*(p/theta1+theta2));

% Test the accurate and approximate criteria
dcrit = double(subs(myFunc,parSym84,par84_long2)); % approximation
% dcrit = double(subs(Bfunction84,parSym84,par84_long2)); % the accurate form

figure(84)
subplot(2,2,1)
d_this = dcrit*1.05;
par84_long2(1) = d_this; 
[~,y] = runModel(@ODE_base_comp,[0.5,0.5],par84_long2,tspan);
title('violates criterion: B/\delta<1')
xlabel('Time (\tau)');ylabel('Density')
axis([0 tspan 0 1])

subplot(2,2,2)
hold on
VectorF(@ODE_base_comp,par84_long2,[0 1 0 1])
ss_val = getSS(ss_comp_84,parSym84,par84_long2); 
for i = 1:length(ss_val)
    if sum(ss_val(i,:)>=0)==2 && sum(abs(imag(ss_val(i,:)))<eps)==2
        plot(real(ss_val(i,1)),real(ss_val(i,2)),'ko')
    end
end
plot(y(:,1),y(:,2))

subplot(2,2,3)
d_this = dcrit*0.95;
par84_long2(1) = d_this;
[~,y] = runModel(@ODE_base_comp,[0.5,0.5],par84_long2,tspan);
title('meets criterion: B/\delta>1')
xlabel('Time (\tau)');ylabel('Density')
axis([0 tspan 0 1])

subplot(2,2,4)
hold on
VectorF(@ODE_base_comp,par84_long2,[0 1 0 1])
par84_long2(1) = d_this;
ss_val = getSS(ss_comp_84,parSym84,par84_long2); 
for i = 1:length(ss_val)
    if sum(ss_val(i,:)>=0)==2 && sum(abs(imag(ss_val(i,:)))<eps)==2
        plot(real(ss_val(i,1)),real(ss_val(i,2)),'ko')
    end
end
plot(y(:,1),y(:,2))

%% Model 85
parSym85 = [d d b b c c 1 1 d0 d0 p];
ss_comp_85 = solve(ODE_base_comp(0,[x1;x2],parSym85)==[0;0],[x1,x2], 'MaxDegree', 3);
par85_long =  [1.2 1.2 100 100 1.1 1.1 1 1 0.1 0.1 1];
par85_long2 = [1.2 1.2 5 5 1.1 1.1 1 1 0.1 0.1 1.1];
%=============B function (approximated)================
Bfunction85 = (b*(1-c*d0)+2)^2/(4*b*c)*p/(p+1);
%======================================================

% verify the criterion
dcrit = double(subs(Bfunction85,parSym85,par85_long2));
figure(85)
subplot(2,2,1)
d_this = dcrit*1.05;
par85_long2(1) = d_this; par85_long2(2) = d_this;
[~,y] = runModel(@ODE_base_comp,[0.5,0.5],par85_long2,tspan);
title('violates criterion: B/\delta<1')
xlabel('Time (\tau)');ylabel('Density')
axis([0 tspan 0 1])

subplot(2,2,2)
hold on
VectorF(@ODE_base_comp,par85_long2,[0 1 0 1])
ss_val = zeros(6,2);
for i = 1:6
    ss_val(i,:) = [double(subs(vpa(ss_comp_85.x1(i)),parSym85,par85_long2)),double(subs(vpa(ss_comp_85.x1(i)),parSym85,par85_long2))];
end
for i = 1:length(ss_val)
    if sum(ss_val(i,:)>=0)==2 && sum(abs(imag(ss_val(i,:)))<eps)==2
        plot(real(ss_val(i,1)),real(ss_val(i,2)),'ko')
    end
end
plot(y(:,1),y(:,2))

subplot(2,2,3)
d_this = dcrit*0.95;
par85_long2(1) = d_this; par85_long2(2) = d_this;
[~,y] = runModel(@ODE_base_comp,[0.5,0.5],par85_long2,tspan);
title('meets criterion: B/\delta>1')
xlabel('Time (\tau)');ylabel('Density')
axis([0 tspan 0 1])

subplot(2,2,4)
hold on
VectorF(@ODE_base_comp,par85_long2,[0 1 0 1])
ss_val = zeros(6,2);
for i = 1:6
    ss_val(i,:) = [double(subs(vpa(ss_comp_85.x1(i)),parSym85,par85_long2)),double(subs(vpa(ss_comp_85.x2(i)),parSym85,par85_long2))];
end
for i = 1:length(ss_val)
    if sum(ss_val(i,:)>=0)==2 && sum(abs(imag(ss_val(i,:)))<eps)==2
        plot(real(ss_val(i,1)),real(ss_val(i,2)),'ko')
    end
end
plot(y(:,1),y(:,2))

%% N-mutualist system
parSym86 = [d,b,c,n]; par86 = [1.5,2,1.1,10];

ss86 = solve((n/c)*x*(1-n*x)-n*d/(b*(n-1)*x+1)*x==0,x);
ss86Sqrt = unique(underSqrt(ss86));
%=============B function ================
eval(['Bfunction86 = solve((',ss86Sqrt{1},')==0,d);'])
%========================================

% verify the criterion
dcrit = double(subs(Bfunction86,parSym86,par86));
figure(86)
subplot(2,1,1)
d_this = dcrit*1.05;
par86(1) = d_this;
runModel(@ODE_nMutualist,ones(10,1)*0.1,par86,tspan);
title('violates criterion: B/\delta<1')
xlabel('Time (\tau)');ylabel('Density')
axis([0 tspan 0 0.2])

subplot(2,1,2)
d_this = dcrit*0.95;
par86(1) = d_this;
runModel(@ODE_nMutualist,ones(10,1)*0.1,par86,tspan);
title('meets criterion: B/\delta>1')
xlabel('Time (\tau)');ylabel('Density')
axis([0 tspan 0 0.2])

