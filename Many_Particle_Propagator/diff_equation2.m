function [Xi_euler, Xi_heun ] = diff_equation2(t, omega, EigVal)
dtau        = t(2) - t(1);
Xi          = zeros(2,2,length(t));
Om          = sqrtm(omega(:,:,1));
Xi(:,:,1)   = Om;
for i = 2:length(t)
    Om2         = omega(:,:,i-1);
    Xi(:,:,i)   = Xi(:,:,i-1) + dtau * (Om2 - Xi(:,:,i-1)*Xi(:,:,i-1));
end
Xi_euler = Xi;
Xi(:,:,1)   = Om;
for i = 2:length(t)
    Om2         = omega(:,:,i-1);
    Om2t        = omega(:,:,i);
    Xi_tmp      = Xi(:,:,i-1) + dtau * (Om2 - Xi(:,:,i-1)*Xi(:,:,i-1));
    Xi(:,:,i)   = Xi(:,:,i-1) + dtau/2 * (( Om2 - Xi(:,:,i-1)*Xi(:,:,i-1)) + (Om2t - Xi_tmp*Xi_tmp));
    
end
Xi_heun = Xi;
end

