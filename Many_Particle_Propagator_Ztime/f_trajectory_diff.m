function diffcurve = f_trajectory_diff(N, traj, z )
    diffcurve = zeros(3, N);
    dz = z(2) - z(1);
    for i = 1:3
        for j = 1:N
            if j == 1
                diffcurve(i, j) = 0;
            elseif j == N
                diffcurve(i, j) = 0;
            else
                diffcurve(i, j) = (traj(i, j + 1) - traj(i, j - 1))/(2 * dz);
            end
        end
    end
    
end