function [f, delta_e, delta_phi, R] = TrajectoryRotationMatrix(e, NoP, Nt)
    %TRAJECTORYROTATIONMATRIX:
    
    delta_e = zeros(NoP, Nt);
    
    for i = 1:Nt-1
        delta_e(:, i) = e(:, i + 1) - e(:, i);
    end

    %Determinng the erpendicular vector
    f               = zeros(NoP, Nt);
    Nominator       = f;
    
    for i = 1:Nt
        n1 = norm(delta_e(:, i));
        n2 = e(:, i)' * delta_e(:, i);
    
        Normalization(i)    = sqrt(abs(n1^2 - n2^2));
        Nominator(:, i)     = delta_e(:, i) - ((e(:, i)' * delta_e(:, i)) * e(:, i));
        if Normalization(i) == 0
            f(:, i) = f(:, i - 1);
        else
            f(:, i) = Nominator(:, i) / Normalization(i);
        end
    end

    % Determining the rotation angle:
    for i = 1:Nt
        if i == 1
            delta_phi(i) = 2 * asin((norm(delta_e(:, i))/2));
        else
            delta_phi(i) = 2 * asin((norm(delta_e(:, i))/2));
        end
    end

    R = zeros(NoP, NoP, Nt);

    for i = 1:Nt
        R1 = eye(NoP) - (e(:,i)*e(:,i)' + f(:,i)*f(:,i)') + (e(:,i)*e(:,i)' + f(:,i)*f(:,i)') * cos((delta_phi(i)));
        R2 = (f(:,i)*e(:,i)' - (e(:,i)*f(:,i)')) * sin((delta_phi(i)));
        R(:,:,i) =  (R1 + R2);
    end

end

