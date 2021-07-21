function [Tauval,Xi_mtx] = diff_equation(omega_t,time,EigVal)

    eigvalues = diag(EigVal);
    EigVal = zeros(2,2);
    for q = 1:2
        EigVal(q,q) = eigvalues(q + 1);
    end
   
    X = zeros(2,2);
    N = length(time);
    for i = 1:N-1
        Tspan = [time(i) time(i+1)];
        omega =  omega_t(:,:,i);
        if i == 1
            X(:,:) = sqrtm(EigVal); 
        else
            X(1,1) = Xi_mtx(1,end);
            X(1,2) = Xi_mtx(2,end);
            X(2,1) = Xi_mtx(3,end);
            X(2,2) = Xi_mtx(4,end);
        end
        
        %[Tau,Xi] = ode113(@(t, xi) [omega(1,1) - (xi(1)); omega(1,2) - (xi(2)); omega(2,1) - (xi(3)); omega(2,2) - (xi(4))], Tspan, [X(1,1) X(1,2) X(2,1) X(2,2)]);
        [Tau,Xi] = ode45(@(t, xi) [omega(1,1) - (xi(1)^2 + xi(2)*xi(3)); omega(1,2) - (xi(1)*xi(2) + xi(2)*xi(4)); omega(1,2) - (xi(3)*xi(1) + xi(3)*xi(4)); omega(2,2) - (xi(4)^2 + xi(3)*xi(2))], Tspan, [X(1,1) X(1,2) X(2,1) X(2,2)]);
        if i == 1
            Xi_mtx = Xi';   %mert elõtte N x 9 es a mátrix volt és én így jobban szeretem
            Tauval = Tau';
        else
            Xi_mtx = [Xi_mtx Xi(2:end-1,:)'];
            Tauval = [Tauval Tau(2:end-1)'];
        end
  
    end
end

