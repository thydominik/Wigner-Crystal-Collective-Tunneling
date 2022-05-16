function [omega] = f_fitting(omega, z_t)
    temp_mtx    = omega(1, 1, 1:50);
    z           = z_t(1:50);
    [gof1, fc1] = f_fitting_omega(z(:), temp_mtx(:));
    omega(1,1,1:25) = omega(1,1,1) + fc1.b .* (z(1:25) + 1) + fc1.c .* (z(1:25) + 1).^2 + fc1.d .* (z(1:25) + 1).^3 + fc1.e .* (z(1:25) + 1).^4;
    
    temp_mtx    = omega(1, 2, 1:50);
    [gof2, fc2] = f_fitting_omega(z(:), temp_mtx(:));
    omega(1,2,1:25) = omega(1,2,1) + fc2.b .* (z(1:25) + 1) + fc2.c .* (z(1:25) + 1).^2 + fc2.d .* (z(1:25) + 1).^3 + fc2.e .* (z(1:25) + 1).^4;
    
    temp_mtx    = omega(2, 1, 1:50);
    [gof3, fc3] = f_fitting_omega(z(:), temp_mtx(:));
    omega(2,1,1:25) = omega(2,1,1) + fc3.b .* (z(1:25) + 1) + fc3.c .* (z(1:25) + 1).^2 + fc3.d .* (z(1:25) + 1).^3 + fc3.e .* (z(1:25) + 1).^4;
    
    temp_mtx    = omega(2, 2, 1:50);
    [gof4, fc4] = f_fitting_omega(z(:), temp_mtx(:));
    omega(2,2,1:25) = omega(2,2,1) + fc4.b .* (z(1:25) + 1) + fc4.c .* (z(1:25) + 1).^2 + fc4.d .* (z(1:25) + 1).^3 + fc4.e .* (z(1:25) + 1).^4;

end

