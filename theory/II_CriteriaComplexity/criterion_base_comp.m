%% This script has
% 1. bifurcation analysis of the base model
% 2. establish the criterion
% 3. establish B/de can predict density and probability


%% find analytical solution to steady states
syms x1 x2 
syms d1 d2 b1 b2 c1 c2 a1 a2 d01 d02 p
syms d b c a d0

ss_comp_1 = solve(ODE_base_comp(0,[x1;x2],[d d b b c c a a 0 0 1])==[0;0],[x1,x2]);

ss_comp_2 = solve(ODE_base_comp(0,[x1;x2],[d1 d2 b b c c 0 0 1])==[0;0],[x1,x2]);
ss_comp_3 = solve(ODE_base_comp(0,[x1;x2],[d d b1 b2 c c 0 0 1])==[0;0],[x1,x2]);
ss_comp_4 = solve(ODE_base_comp(0,[x1;x2],[d d b b c1 c2 0 0 1])==[0;0],[x1,x2]);
ss_comp_5 = solve(ODE_base_comp(0,[x1;x2],[d1 d2 b1 b2 c1 c2 0 0 p])==[0;0],[x1,x2]);

ss_comp_6 = solve(ODE_base_comp(0,[x1;x2],[d d b b c c d0 d0 1])==[0;0],[x1,x2]);
ss_comp_7 = solve(ODE_base_comp(0,[x1;x2],[d1 d2 b b c c d0 d0 1])==[0;0],[x1,x2]);
bene7 = solve((4*b + b^2 - 4*b*c*d0 - 4*b*c*d1 - 4*b*c*d2 + b^2*c^2*d0^2 - 2*b^2*c*d0 + 4)==0,c*(d1+d2));

ss_comp_8 = solve(ODE_base_comp(0,[x1;x2],[d1 d2 b1 b2 c c d0 d0 1])==[0;0],[x1,x2]);

