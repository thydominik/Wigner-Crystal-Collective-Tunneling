function [diff_trajectory] = traj_diff(particle_number, N, lim, cf, eqp)
    diff_trajectory = zeros(particle_number, N);   
    %analitikus megoldás:
    x = linspace(-lim, lim, N);
   
    diff_trajectory(1,:) = cf(2,1)*cf(3,1)./(cosh(cf(4,1) + cf(3,1).*x).^2);
    diff_trajectory(2,:) = cf(2,2)*cf(3,2)./(cosh(cf(4,2) + cf(3,2).*x).^2);
    diff_trajectory(3,:) = cf(2,3)*cf(3,3)./(cosh(cf(4,3) + cf(3,3).*x).^2);
        
end

