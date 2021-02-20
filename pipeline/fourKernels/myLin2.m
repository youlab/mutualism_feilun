function G = myLin2(U,V)

global k1 k2

G = U(:,1:2)*V(:,1:2)'+ k1 + k2*U(:,3)*V(:,3)'; 

end