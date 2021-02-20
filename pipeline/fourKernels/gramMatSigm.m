% gram matrix of sigmoidal kernel
% setup the data vector

[X,Y]=meshgrid(linspace(-10,10,21),linspace(-10,10,21));

%% 
global k1 k2
k1 = -3;
k2 = 3;

%
% figure(1)
% imagesc(myLin2(X,Y))
% colorbar
% %%
% figure(2)
% imagesc(mySigm2(X,Y))
% colorbar

%%
this_point = [2,2];
linK = myCub2V([X(:),Y(:)],this_point);
figure(3)
imagesc(reshape(linK,21,21))
colorbar
set(gca,'xtick',1:21,'xticklabel',-10:10)
set(gca,'ytick',1:21,'yticklabel',-10:10)

%%
this_point = [2,2];
% linK = mySigm2V([X(:),Y(:)],this_point);
linK = mySigm2V(this_point,this_point);
figure(4)
imagesc(reshape(linK,21,21))
colorbar
set(gca,'xtick',1:21,'xticklabel',-10:10)
set(gca,'ytick',1:21,'yticklabel',-10:10)

% diag(([X(:),Y(:)]-this_point)*([X(:),Y(:)]-this_point)')