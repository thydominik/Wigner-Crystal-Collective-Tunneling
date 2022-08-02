function [C] = f_curvature(e, tauspace, z_time, v)
    C = zeros(3,length(tauspace(1,1,:)));

    %numerically derivating the tauspace vectors
    for i = 2: length(C)-1
        %dt here is 2 times dz
        dt              = z_time(3) - z_time(1);
        derivative1     = (tauspace(:,1,i+1) - tauspace(:,1,i-1))/(dt); %- 2*tauspace(:,1,i))/(2*dt);
        derivative2     = (tauspace(:,2,i+1) - tauspace(:,2,i-1))/(dt); %- 2*tauspace(:,2,i))/(2*dt);
        derivative3     = (tauspace(:,3,i+1) - tauspace(:,3,i-1))/(dt); %- 2*tauspace(:,3,i))/(2*dt);
    
        der = [derivative1 derivative2 derivative3];
        
        C(1,i) = (der(:,1)' * tauspace(:,1,i)); %e(:,i);  %-vel
        C(2,i) = (der(:,2)' * tauspace(:,2,i));
        C(3,i) = (der(:,3)' * tauspace(:,3,i));
        
    end
end
