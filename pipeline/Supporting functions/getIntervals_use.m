function [Variances_,B_u_,B_l_,B_m_,percInt_]=getIntervals_use(X,Y,k1_,k2_,ker,grid,BoxC,n)

% This function runs the simple boostrapping function (conInt_Boot_use()) 
% with a set of kernel parameters. 

global k1 k2
Variances_ = cell(length(k1_),length(k2_));
B_u_ = cell(length(k1_),length(k2_));
B_l_ = cell(length(k1_),length(k2_));
B_m_ = cell(length(k1_),length(k2_));
percInt_ = cell(length(k1_),length(k2_));

for i = 1:length(k1_)
    for j = 1:length(k2_)
        k1 = k1_(i);
        k2 = k2_(j);
        eval(['[B_u,B_l,B_m,Variances]=conInt_Boot_use(X,Y,n,grid,ker,@',ker,'V,BoxC);'])
        
        Variances_{i,j}= Variances;
        B_u_{i,j} = B_u;
        B_l_{i,j} = B_l;
        B_m_{i,j} = B_m;
        percInt_{i,j}= abs(B_u-B_l)./abs(B_m);
    end
end
