function [ Xi_heun, temppp ] = f_diff_equation(t, omega, r, NoP)
    
    
    dtau        = t(2) - t(1);
    Xi          = zeros(NoP - 1, NoP - 1,length(t));
    Om          = sqrtm(omega(:,:,1));
    Xi(:,:,1)   = Om;
    % for i = 2:length(t)-1
    %     dtau = atanh(t(i)) - atanh(t(i-1));
    %     Om2         = omega(:,:,i-1);
    %     Xi(:,:,i)   = Xi(:,:,i-1) + dtau * (Om2 - Xi(:,:,i-1)*Xi(:,:,i-1));
    % end
    
    Xi_euler = Xi;
    Xi(:,:,1)   = Om;
    Xi(:,:,2)   = sqrtm(omega(:, :, 2));
    for i = 3:length(t)-1
    
        
        dtau = t(3) - t(2);
        Om2         = sqrtm(omega(:,:,i-1));
        Om2         = Om2 * Om2;
        Om2t        = omega(:,:,i);
        
        if (Om2 - Xi(:,:,i-1)*Xi(:,:,i-1)) == zeros(NoP-1, NoP-1)
            Xi_tmp = Xi(:,:,i-1);
            disp('okay')
        else
            Xi_tmp      = Xi(:,:,i-1) + dtau * r/(1 - t(i-1)^2) * (Om2 - Xi(:,:,i-1)*Xi(:,:,i-1));
        end
        
        if ( Om2 - Xi(:,:,i-1)*Xi(:,:,i-1)) ==  zeros(NoP-1, NoP-1) 
            Xi_i_1 = zeros(NoP-1, NoP-1);
            disp('okay')
        else
            Xi_i_1 = dtau/2 * (r/(1 - t(i-1)^2) *( Om2 - Xi(:,:,i-1)*Xi(:,:,i-1)));
        end
        
        if (Om2t - Xi_tmp*Xi_tmp) ==  zeros(NoP-1, NoP-1)
            Xi_i_2 = zeros(NoP-1, NoP-1);
        else
            Xi_i_2 = dtau/2 * r/(1 - t(i)^2) *(Om2t - Xi_tmp*Xi_tmp);
        end
        
        Xi(:,:,i)   = Xi(:,:,i-1) + Xi_i_1 + Xi_i_2;
        
        temppp(: ,:,i)      = sqrtm(Om2) - Xi(:, :, i-1);
    end
    %Xi(:,:,end) = sqrtm(om(:,:,end));
    Xi_heun = Xi;
end


