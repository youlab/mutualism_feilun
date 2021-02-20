function dydt = ODE_base_comp(time,y,params)

    d1=params(1);
    d2=params(2);
    b1=params(3);
    b2=params(4);
    c1=params(5);
    c2=params(6);
    a1=params(7);
    a2=params(8);    
    d01=params(9);
    d02=params(10);
    p=params(11);

    dydt(1:2,:,:) = [  (1/c1)*y(1,:,:).*(1-y(1,:,:)-a2*y(2,:,:))-d1./(b2*y(2,:,:)+1).*y(1,:,:)-d01.*y(1,:,:);
                     p*(1/c2)*y(2,:,:).*(1-y(2,:,:)-a1*y(1,:,:))-d2./(b1*y(1,:,:)+1).*y(2,:,:)-d02.*y(2,:,:);
            ];
    

