function [S_int] = f_Action(q0, r, a, div, z)
    dz = z(2) - z(1);
    
    der = zeros(1,div);
    der(1) = (q0(2)-q0(1))/dz;

    for i = 2:div-1
        %der(i) = ((q0(i+1) - q0(i-1))/(2*dz));
        der(i) = ((q0(i+1) - q0(i))/(dz));
    end
    der(div) = ((q0(div)-q0(div-1))/(dz));
    
    func = zeros(1,div);
    for k=1:div
        pre = (1 - (z(k)^2))/r;
        func(k) = pre*(1/2) * der(k)^2 + (1/pre)*(1/4)*(q0(k)^2 - a)^2;
    end
    
    %trapezoid rule
    actionn = 0;
    for L = 2:div
        actionn = actionn + (dz * (func(L-1)+func(L))/2);   
    end
    S_int = actionn;
end