function [Biass_,Variances_,Biass2_,Variances2_,RSD_,STD_,Amp_]=getIntervals(X,Y,k1_,k2_,ker,grid,BoxC,n,B_real)

% This function runs the boostrapping function (conInt_Boot()) with a set
% of kernel parameters. 

global k1 k2
Variances_ = cell(length(k1_),length(k2_));
Biass_ = cell(length(k1_),length(k2_));
Variances2_ = cell(length(k1_),length(k2_));
Biass2_ = cell(length(k1_),length(k2_));
RSD_ = cell(length(k1_),length(k2_));
STD_ = cell(length(k1_),length(k2_));
Amp_ = cell(length(k1_),length(k2_));

for i = 1:length(k1_)
    for j = 1:length(k2_)
        k1 = k1_(i);
        k2 = k2_(j);
        eval(['[Biass,Variances,Biass2,Variances2,RSD,STD,Amp]=conInt_Boot(X,Y,n,grid,ker,@',ker,'V,BoxC,B_real);'])
        
        Biass_{i,j} = Biass;
        Variances_{i,j}= Variances;
        Biass2_{i,j} = Biass2;
        Variances2_{i,j}= Variances2;   
        RSD_{i,j} = RSD;
        STD_{i,j} = STD;
        Amp_{i,j} = Amp;
    end
end
