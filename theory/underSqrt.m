function outputs = underSqrt(ss)

% This function returns the expressions under the square root of the steady
% states.

% INPUTS
% ss: steady state solutions of one population

n = length(ss);
ind = 1;
outputs = cell(1,1);
for i = 1:n
    sol = char(ss(i));

    a = strfind(sol,')^(1/2)'); % need to be refined further
    a = a+1;

    if ~isempty(a)&&length(a)==1
        lp = find(sol=='(');
        rp = find(sol==')');
        flag = 1;
        k = a-1; % starting index of what is before '^(1/2)'
        while flag ~=0 && k>0
            k = k-1; 
            if ismember(k,lp)
                flag = flag - 1;
            elseif ismember(k,rp)
                flag = flag + 1;
            end   
        end
        outputs{ind} = sol(k:(a-1));
        
        ind = ind + 1;
    end
end

end