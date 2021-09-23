function [BB] = f_B_matrix(q, a, tauspace)
    %first just the potential without the interaction:
    eta = 18.813;

    %creating the second derivative diagonal potential elements for every point
    %in the trajectory
    xxV =  (a + 3*(q(1,:).^2));
    yyV =  (a + 3*(q(2,:).^2));
    zzV =  (a + 3*(q(3,:).^2));

    %intercation terms (diag)
    xxU = eta*2*(1./((q(2,:) - q(1,:)).^3) + 1./((q(3,:) - q(1,:)).^3));
    yyU = eta*2*(1./((q(3,:) - q(2,:)).^3) + 1./((q(2,:) - q(1,:)).^3));
    zzU = eta*2*(1./((q(3,:) - q(1,:)).^3) + 1./((q(3,:) - q(2,:)).^3));

    %interaction terms (off-diag)
    xyU = -eta*(2./((q(2,:) - q(1,:)).^3));
    xzU = -eta*(2./((q(3,:) - q(1,:)).^3));
    yzU = -eta*(2./((q(3,:) - q(2,:)).^3));

    V = zeros(3,3,length(q));
    V(1,1,:) = xxV + xxU;
    V(2,2,:) = yyV + yyU;
    V(3,3,:) = zzV + zzU;

    V(1,2,:) = xyU;
    V(2,1,:) = V(1,2,:);
    V(1,3,:) = xzU;
    V(3,1,:) = V(1,3,:);
    V(2,3,:) = yzU;
    V(3,2,:) = V(2,3,:);

    %this is 3 by 3 by N_division size tensor
    B =  V;

    BB = zeros(2,2,length(B(1,1,:)));
    for i = 1:2
        for j = 1:2
            for k = 1:length(tauspace(1,1,:))
                temp        =  B(:,:,k) * tauspace(:,j+1,k);
                %temp        =  tauspace(:,i+1,k)' * B(:,:,k);
                temp2       =  tauspace(:,i+1,k)' * temp;
                %temp2       =  temp * tauspace(:,j+1,k);
                BB(i,j,k)   =  temp2;
            end
        end
    end
end

