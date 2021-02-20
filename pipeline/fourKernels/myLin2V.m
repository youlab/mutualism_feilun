function G = myLin2V(U,V)

global k1

G = U(:,1:2)*V(:,1:2)'+ k1 ; 

end