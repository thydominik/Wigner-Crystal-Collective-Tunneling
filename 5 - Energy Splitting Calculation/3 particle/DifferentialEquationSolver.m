function [Xi] = DifferentialEquationSolver(z, dz, NoP, Nt, Omega, r)
%DIFFERENTIALEQUATIONSOLVER:
Xi = zeros(NoP - 1, NoP - 1, Nt);

% Initial condition:
Xi(:, :, 1) = sqrtm(Omega(:, :, 1));
Xi(:, :, 2) = sqrtm(Omega(:, :, 2));
Xi(:, :, 3) = sqrtm(Omega(:, :, 3));
Xi(:, :, 4) = sqrtm(Omega(:, :, 4));
Xi(:, :, 5) = sqrtm(Omega(:, :, 5));
Xi(:, :, 6) = sqrtm(Omega(:, :, 6));
% Heun's method:
% for time_ind = 3:Nt
%     
%     Omega2  = Omega(:, :, time_ind - 1);
%     Omega2t = Omega(:, :, time_ind);
% 
%     if (Omega2 - Xi(:, :, time_ind - 1) * Xi(:, :, time_ind - 1)) == zeros(NoP - 1, NoP - 1)
%         Xi_temp1 = Xi(:, :, time_ind - 1);
%         disp('ok')
%     else
%         Xi_temp1 = Xi(:, :, time_ind - 1) + 0.5 * dz * r/(1 - z(time_ind - 1)^2) * (Omega2 - Xi(:, :, time_ind - 1)*Xi(:, :, time_ind));
%     end
% 
%     if (Omega2t - Xi_temp1*Xi_temp1) == zeros(NoP-1, NoP-1)
%         Xi_temp2 = zeros(NoP-1, NoP-1);
%     else
%         Xi_temp2 = 0.5 * dz * r/(1 - z(time_ind)^2) * (Omega2t - Xi_temp1*Xi_temp1);
%     end
%     
%     Xi(:, :, time_ind) = Xi(:, :, time_ind - 1) + 0.5 * Xi_temp2 + 0.5* (Xi_temp1);
% end

for time_ind = 7:Nt
    F1 = dz * r /(1 - z(time_ind - 1)^2) ;
    F2 = Omega(:, :, time_ind - 1) - (Xi(:, :, time_ind - 1) * Xi(:, :, time_ind - 1));
    Xi_tilde = Xi(:, :, time_ind - 1) + F1 * F2;


    F1tm = dz * r /(1 - z(time_ind - 1)^2);
    F1tp = dz * r /(1 - z(time_ind)^2);
    F2tm = (Omega(:, :, time_ind - 1) - (Xi(:, :, time_ind - 1) * Xi(:, :, time_ind - 1)));
    F2tp = (Omega(:, :, time_ind) - (Xi_tilde * Xi_tilde));
    Xi(:, :, time_ind) = Xi(:, :, time_ind - 1) + 0.5 * (F1tm * F2tm + F1tp * F2tp);
    

end

% for time_ind = 3:2:Nt
%     A = Omega(:, :, time_ind - 1) - Xi(:, :, time_ind)*Xi(:, :, time_ind);
%     B = Omega(:, :, time_ind) - A*A;
%     C = Omega(:, :, time_ind) - B*B;
%     D = Omega(:, :, time_ind + 1) - 
% end


end

