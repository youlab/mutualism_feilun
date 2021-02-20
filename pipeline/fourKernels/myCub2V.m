function G = myCub2V(U,V)

global k1
G = (U(:,1:2)*V(:,1:2)'+ k1).^3 ; 

% G = (U(:,2)*V(:,2)'+ k1).^2 + k2*(U(:,3)-V(:,3)).^2'; 

end