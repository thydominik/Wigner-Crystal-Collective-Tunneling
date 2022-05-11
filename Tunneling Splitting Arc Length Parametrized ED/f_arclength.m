function [chiS, LS, VS] = f_arclength(traj, a, eta, zz)
%F_ARCLENGTH
N = length(traj(:, 1));
LS(1) = 0;
for i = 2:length(traj)
    temp = 0;
    for p = 1:N
        temp = temp + (traj(p, i) - traj(p, i - 1))^2;
    end
    S(i) = sqrt(temp);
    LS(i) = LS(i-1) + S(i);
end

for i = 1:length(traj)
    temp = 0;
    for p = 1:N
        temp = temp + 0.25 * (traj(p, i)^2 + a)^2;
        for q = (p+1):N
            temp = temp + eta/abs(traj(p, i) - traj(q, i));
        end
    end
    VS(i) = temp;
end
chiS = 0;
end

