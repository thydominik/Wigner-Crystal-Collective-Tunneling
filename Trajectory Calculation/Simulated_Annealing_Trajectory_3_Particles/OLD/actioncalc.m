function action = actioncalc(pos,r,a,rs,N,z,dz,shift)
    action = 0;
    %I denoted the interaction part with Q
    Q = zeros(N,1);
    for i = 1:N
        pre = (1 - z(i)^2)/r;
        Q(i) = (rs/pre) * abs((1/(pos(2,i) - pos(1,i))) + (1/(pos(3,i) - pos(1,i))) + (1/(pos(3,i) - pos(2,i)))); 
    end
    
    F1 = zeros(N,1);
    F2 = zeros(N,1);
    F3 = zeros(N,1);
    der1 = zeros(N,1);
    der2 = zeros(N,1);
    der3 = zeros(N,1);
    
    der1(1) = (pos(1,2) - pos(1,1))/dz;
    der2(1) = (pos(2,2) - pos(2,1))/dz;
    der3(1) = (pos(3,2) - pos(3,1))/dz;
    
    for i = 2:N-1
        der1(i) = (pos(1,i) - pos(1,i-1))/(dz);
        der2(i) = (pos(2,i) - pos(2,i-1))/(dz);
        der3(i) = (pos(3,i) - pos(3,i-1))/(dz);
    end
    
    der1(N) = (pos(1,N) - pos(1,N - 1))/dz;
    der2(N) = (pos(2,N) - pos(2,N - 1))/dz;
    der3(N) = (pos(3,N) - pos(3,N - 1))/dz;
    
    for i = 1:N
        pre = (1 - z(i)^2)/r;
        
        F1(i) = pre/2 * der1(i)^2 + (1/pre) * 0.25 * (pos(1,i)^2 + a)^2;
        F2(i) = pre/2 * der2(i)^2 + (1/pre) * 0.25 * (pos(2,i)^2 + a)^2;
        F3(i) = pre/2 * der3(i)^2 + (1/pre) * 0.25 * (pos(3,i)^2 + a)^2;
    end
           
    %trapezoid rule
    for L = 2:N
       action = action + (dz/2 * (F1(L-1) + F1(L) + F2(L-1) + F2(L) + F3(L-1) + F3(L) + Q(L-1) + Q(L))); 
    end
    action = action;
end