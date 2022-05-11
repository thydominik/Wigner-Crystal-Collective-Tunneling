function [prop,Trace] = prefactor( Xi_matrix, EigValu, tau_time, y_time, a )
    %first the factor before the exponential
   eigvalues = diag(EigValu(:,:,1));
   EigVal = zeros(2,2);
   for q = 1:2
       EigVal(q,q) = eigvalues(q + 1);
   end
%    EigVal = (EigValu(:,:,1));
   omega_determinant    = det(sqrtm(EigVal));
   disp(['omega_0 determinant: ', num2str(omega_determinant)])
   
   xi_determinant       = det(abs(Xi_matrix(round(length(y_time)/2))));
   disp(['xi(0) determinant: ', num2str(xi_determinant)])
   
   prefactor            = sqrt((omega_determinant / xi_determinant));
   disp(['prefactor: ', num2str(prefactor)])
   
   % now the exponential and integration
   
   Trace = zeros(1,length(tau_time) ); %-1 becouse Tau index is the first non negative time marker
   for i = 1:length(tau_time)/2
       Trace(i) = trace(sqrt(EigVal) - Xi_matrix(:,:,i));
   end
   
   figure(20)
   clf(figure(20))
   hold on
   title('Trace from T_0 to 0')
   plot(tau_time, smooth(Trace))
   plot(tau_time, zeros(1,length(tau_time)))
   %scatter(y_time, zeros(1,length(y_time)))
   xlim([tau_time(1) 0])
   hold off
   
   %integration (midpoint)
   Integ = 0;
   for i = 1:length(tau_time)/2
       dt = tau_time(i+1) - tau_time(i);
       Integ = Integ + Trace(i) * dt;
   end
   disp(['Integral: ', num2str(Integ)])
   prop = prefactor * exp(Integ);
   disp(['prefactor: sqrt term x exp= ', num2str(prop)])
   
   
   
   Action = 1.136*(-a-4.33333)^(1.483); %0.924*a^(1.49999);
   omega = sqrt(EigValu(1,1));
   splitt = sqrt(2)^2 * omega * sqrt(abs(a)) * sqrt(omega/pi) * exp(-Action);
   prop = prop * splitt;
   disp(['prefactor: 1D * sqrt term x exp= ', num2str(prop)])
end

