function [Xi_euler, Xi_heun ] = f_diff_equation(t, omega)
omega(1,1,:) = smooth(omega(1,1,:));
omega(1,2,:) = smooth(omega(1,2,:));
omega(2,1,:) = smooth(omega(2,1,:));
omega(2,2,:) = smooth(omega(2,2,:));

dtau        = t(2) - t(1);
Xi          = zeros(2,2,length(t));
Om          = sqrtm(omega(:,:,1));
Xi(:,:,1)   = Om;
% for i = 2:length(t)-1
%     dtau = atanh(t(i)) - atanh(t(i-1));
%     Om2         = omega(:,:,i-1);
%     Xi(:,:,i)   = Xi(:,:,i-1) + dtau * (Om2 - Xi(:,:,i-1)*Xi(:,:,i-1));
% end
Xi_euler = Xi;
Xi(:,:,1)   = Om;
for i = 2:length(t)-1
    if t(i-1) == -1
        t(i - 1) = -1 + 10^-15;
    end
    
    dtau = t(i) - t(i-1);
    Om2         = omega(:,:,i-1);
    Om2t        = omega(:,:,i);
    Xi_tmp      = Xi(:,:,i-1) + dtau * 1/(1 - t(i-1)^2) * (Om2 - Xi(:,:,i-1)*Xi(:,:,i-1));
    Xi(:,:,i)   = Xi(:,:,i-1) + dtau/2 * (1/(1 - t(i-1)^2) *( Om2 - Xi(:,:,i-1)*Xi(:,:,i-1)) + 1/(1 - t(i)^2) *(Om2t - Xi_tmp*Xi_tmp));
    
end
Xi_heun = Xi;
end


