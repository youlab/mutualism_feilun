function ss=getSS(ssSolve,symbols,params)

% This function outputs numerical values of steady statesf with a certain
% set of parameters. 

% INPUTS
% ssSolve: symbolic expressions of steady states
% symbols: a vector of all symbols in the expressions
% params: numerical values of parameters

n = length(ssSolve.x1);
ss = zeros(n,2);

for i = 1:n
    ss(i,1) = double(subs(ssSolve.x1(i),symbols,params));
    ss(i,2) = double(subs(ssSolve.x2(i),symbols,params));
end

end