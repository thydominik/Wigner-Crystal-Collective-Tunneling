function  R = f_rotation_matrix(N, e, f, delta_phi, NoP)
    R = zeros(NoP,NoP,N);
    for i = 1:N
        R1 = eye(5) - (e(:,i)*e(:,i)' + f(:,i)*f(:,i)') + (e(:,i)*e(:,i)' + f(:,i)*f(:,i)') * cos((delta_phi(i)));
        R2 = (f(:,i)*e(:,i)' - (e(:,i)*f(:,i)')) * sin((delta_phi(i)));
        R(:,:,i) =  (R1 + R2);
    end
end

