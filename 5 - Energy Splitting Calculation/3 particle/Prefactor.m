function [prop, Trace1, Trace2, Integ1, Integ2, Integ3] = Prefactor(Xi, Omega, EigVals, z, dz, Alpha, R, NoP, Nt)
    %PREFACTOR:
    EigVals = EigVals(2:end, 2:end);
    
    OmegaDeterminant = det(sqrtm(EigVals));
    disp(['omega_0 determinant: ', num2str(OmegaDeterminant)])
    
    XiDeterminant    = det(abs(Xi(:,:,round(Nt/2))));
    disp(['xi(0) determinant: ', num2str(XiDeterminant)])
    
    FractionPrefactor    = sqrt((OmegaDeterminant / XiDeterminant));
    disp(['prefactor: ', num2str(FractionPrefactor)])
    
    Trace1 = zeros(1, Nt);
    %calculate trace way 1:
    for i = 1:Nt/2
        if i == 1
            Trace1(i) = 0;
        elseif i == 2
            Trace1(i) = trace(Xi(:, :, 1) - Xi(:, :, i))/1;
        else
            Trace1(i) = trace(Xi(:, :, 1) - Xi(:, :, i));
        end
    end
    
    %Calculate trace way 2:
    for i = 1:Nt/2
        if i == 1
            Trace2(i) = 0;
        else
            Trace2(i) = trace(sqrtm(EigVals) - Xi(:, :, i));
        end
    end
    

    % integration:
    %integration (midpoint)
   Integ1 = 0;
   for i = 2:Nt/2
       dz       = z(i) - z(i-1);
       Integ1   = Integ1 + (Trace1(i) + Trace1(i - 1)) * dz/2 * R/(1 - z(i)^2);
   end
   %Integ1 = Integ1 + (-z(i)) * Trace1(i);
   
   Integ2 = 0;
   for i = 2:Nt/2
       dz       = z(i+1) - z(i);
       Integ2   = Integ2 + (Trace2(i) + Trace2(i - 1)) * dz/2 * R/(1 - z(i)^2);
   end
   Integ2 = Integ2 + (-z(i)) * Trace2(i);
   
   Integ3 = 0;
   Trace3 = Trace1;
   Trace3(10:100) = (Trace1(10:100) + Trace2(10:100))/2;
   for i = 2:Nt/2
       dz = z(i+1) - z(i);
       Integ3 = Integ3 + (Trace3(i) + Trace3(i - 1)) * dz/2 * R/(1 - z(i)^2);
   end
   Integ3 = Integ3 + (-z(i)) * Trace3(i);
   disp(['Integral 1: ', num2str(Integ1)])
   disp(['Integral 2: ', num2str(Integ2)])
   disp(['Integral 3: ', num2str(Integ3)])
   prop = FractionPrefactor * exp(Integ1); %előbb Integ3 volt beírva
   disp(['prefactor: sqrt term x exp= ', num2str(prop)])
end

