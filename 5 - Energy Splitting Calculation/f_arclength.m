function [chiS, LS, VS] = f_arclength(traj, a, eta, zz)
%F_ARCLENGTH
a = -a;
LS(1) = 0;
for i = 2:length(traj)
    x = traj(1, i) - traj(1, i-1);
    y = traj(2, i) - traj(2, i-1);
    z = traj(3, i) - traj(3, i-1);
    S(i) = sqrt(x^2 + y^2 + z^2);
    LS(i) = LS(i-1) + S(i);
end

for i = 1:length(traj)
    VS(i) = 0.25 * ((traj(1, i)^2 + a)^2 + (traj(2, i)^2 + a)^2 + (traj(3, i)^2 + a)^2);
    VS(i) = VS(i) + eta * (1/abs(traj(1, i) - traj(2, i)) + 1/abs(traj(1, i) - traj(3, i)) + 1/abs(traj(2, i) - traj(3, i)));
end
chiS = 0;
end

