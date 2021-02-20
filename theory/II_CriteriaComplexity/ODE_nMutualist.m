function dydt = ODE_nMutualist(time,y,params)

    d = params(1);
    b = params(2);
    c = params(3);
    n = params(4);

    ysum = sum(y);
    
    b_matrix = ones(n,n)*b;
    b_matrix(eye(n)==1) = 0;
    
    Binter = b_matrix*y;
    
    dydt = 1/c*(1-ysum).*y-d./(Binter+1).*y;
        