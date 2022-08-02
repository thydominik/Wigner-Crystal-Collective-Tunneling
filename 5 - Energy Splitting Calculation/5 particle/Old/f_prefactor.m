function [prop, Trace1, Trace2] = f_prefactor( Xi_matrix, EigValu, tau_time, z_time, a, r, NoP)
   %first the factor before the exponential

   eigvalues = diag(EigValu(:,:,1));
   EigVal = zeros(NoP - 1, NoP - 1);

   for q = 1:NoP-1
       EigVal(q,q) = eigvalues(q + 1);
   end
%    EigVal = (EigValu(:,:,1));
   omega_determinant    = det(sqrtm(EigVal));
   disp(['omega_0 determinant: ', num2str(omega_determinant)])
   
   xi_determinant       = det(abs(Xi_matrix(:, :, round(length(z_time)/2))));
   disp(['xi(0) determinant: ', num2str(xi_determinant)])
   
   prefactor            = sqrt((omega_determinant / xi_determinant));
   disp(['prefactor: ', num2str(prefactor)])
   
   % now the exponential and integration
   
   Xi_matrix(1:(NoP-1), 1:(NoP -1), end-1) = sqrtm(EigVal);
   
   
   Trace = zeros(1,length(tau_time) ); %-1 becouse Tau index is the first non negative time marker
   for i = 1:length(tau_time)/2
       if i == 1
            Trace1(i) = 0;
       elseif i == 2
            Trace1(i) = trace(sqrt(EigVal) - Xi_matrix(:,:,end - (i+1) + 1))/2;
       else
            Trace1(i) = trace(sqrt(EigVal) - Xi_matrix(:,:,end - i + 1));
       end
   end
   for i = 1:length(tau_time)/2
       if i == 1
            Trace2(i) = 0;
       else
            Trace2(i) = trace(sqrt(EigVal) - Xi_matrix(:,:,i));
       end
   end
 
   %integration (midpoint)
   Integ1 = 0;
   for i = 2:length(tau_time)/2
       dt       = tau_time(i+1) - tau_time(i);
       Integ1   = Integ1 + (Trace1(i) + Trace1(i - 1)) * dt/2 * r/(1 - z_time(i)^2);
   end
   Integ1 = Integ1 + (-tau_time(i)) * Trace1(i);
   
   Integ2 = 0;
   for i = 2:length(tau_time)/2
       dt       = tau_time(i+1) - tau_time(i);
       Integ2   = Integ2 + (Trace2(i) + Trace2(i - 1)) * dt/2 * r/(1 - z_time(i)^2);
   end
   Integ2 = Integ2 + (-tau_time(i)) * Trace2(i);
   
   Integ3 = 0;
   Trace3 = Trace1;
   Trace3(10:100) = (Trace1(10:100) + Trace2(10:100))/2;
   for i = 2:length(tau_time)/2
       dt = tau_time(i+1) - tau_time(i);
       Integ3 = Integ3 + (Trace3(i) + Trace3(i - 1)) * dt/2 * r/(1 - z_time(i)^2);
   end
   Integ3 = Integ3 + (-tau_time(i)) * Trace3(i);
   disp(['Integral 1: ', num2str(Integ1)])
   disp(['Integral 2: ', num2str(Integ2)])
   disp(['Integral 3: ', num2str(Integ3)])
   prop = prefactor * exp(Integ1); %előbb Integ3 volt beírva
   disp(['prefactor: sqrt term x exp= ', num2str(prop)])
end

