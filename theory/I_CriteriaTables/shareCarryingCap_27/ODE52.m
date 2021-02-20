function dydt = ODE52(time,y,params)

    d=params(1);
    b=params(2);
    c=params(3);
    b2=params(4);
    
    dydt = [(b2*y(2,:,:)./(y(2,:,:)+1/b)).*y(1,:,:).*(1-y(1,:,:)-y(2,:,:))-d*c.*y(1,:,:);
            (b2*y(1,:,:)./(y(1,:,:)+1/b)).*y(2,:,:).*(1-y(1,:,:)-y(2,:,:))-d*c.*y(2,:,:);
            ];
end