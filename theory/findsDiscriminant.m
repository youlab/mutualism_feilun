function discrim = findsDiscriminant(ss)

% This function returns the expressions under the square root of the steady
% states.

% INPUTS
% ss: steady state solutions of one population
syms z
n = length(ss);

for i = 1:n
    
    sol = char(ss(i));
    l_sol = length(sol);
    rootP = strfind(sol,'root');
    
    inRoot = sol((rootP+5):(l_sol-7));
     
    inRootSym = sym(inRoot);
    coeffPoly = coeffs(inRootSym,z);
    a4 = coeffPoly(1);a3 = coeffPoly(2);a2 = coeffPoly(3);a1 = coeffPoly(4);

    discrim = simplify(18*a1*a2*a3*a4-4*a2^3*a4+a2^2*a3^2-4*a1*a3^3-27*a1^2*a4^2);

end

end