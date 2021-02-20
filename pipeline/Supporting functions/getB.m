function [B] = getB(X,Y,svmModel,grids,ker,Xtot,Ytot)

% This function calculate the B landscape using an input data and
% extrapolate on a meshgrid. 

% if Xtot and Ytot are included in the function input, grids have to be the
% same as Xtot(1:2,:).

% X, Y are the input data used for training svmModel (see below)
% svmModel is the trained classification SVM model using ker (see below)
% grid is the meshgrid the calculated B will be on
% ker is a specific set of kernel parameter
% Xtot,Ytot are the total known data which are used to determine whether a
%   specific B needs to be flipped

global k2

n = length(grids);

isSV = svmModel.IsSupportVector;
alph = svmModel.Alpha;
biasSVM = svmModel.Bias;

B = zeros(n,1);

for i = 1:n
    B(i)= -(sum(alph.*ker(X(isSV,1:end-1),grids(i,:)).*Y(isSV))+biasSVM)./(k2*(sum(alph.*X(isSV,end).*Y(isSV))));
end

if nargin>5
    % this function checks the directionality of the calibrated B. If it is in
    % opposite direction, flip it to the correct dirrection.

    if sum(sign(B-Xtot(:,end)).*Ytot)<0
        B = 2*Xtot(:,end)-(B);
    end
end

end