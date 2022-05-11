function [tauspace] = Tau_space(R,tau)
    tauspace = zeros(3,3,length(R(3,3,:)));
    temp_tau1 = tau(:,1);
    temp_tau2 = tau(:,2);
    temp_tau3 = tau(:,3);
    
    tauspace(:,1,1)     = temp_tau1;
    tauspace(:,2,1)     = temp_tau2;
    tauspace(:,3,1)     = temp_tau3;
    

    for i = 2:length(R(3,3,:))
        temp_tau1           = R(:,:,i)*temp_tau1;
        temp_tau2           = R(:,:,i)*temp_tau2;
        temp_tau3           = R(:,:,i)*temp_tau3;
        tauspace(:,1,i)     = temp_tau1;
        tauspace(:,2,i)     = temp_tau2;
        tauspace(:,3,i)     = temp_tau3;
%         figure(101)
%         hold on
%         quiver3(0,0,0,tauspace(1,1,i),tauspace(2,1,i),tauspace(3,1,i))
%         hold off       
%         pause
    end
    
end

