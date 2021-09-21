function  R = f_rotation_matrix(N, e, f, delta_phi)
    R = zeros(3,3,N);
    for i = 1:N
        R1 = eye(3) - (e(:,i)*e(:,i)' + f(:,i)*f(:,i)') + (e(:,i)*e(:,i)' + f(:,i)*f(:,i)') * cos((delta_phi(i)));
        R2 = (e(:,i)*f(:,i)' - (f(:,i)*e(:,i)')) * sin((delta_phi(i)));
        R(:,:,i) =  (R1 + R2)';
    end
end

