function [chiS, LS, VS] = f_arclength(traj, a, eta, zz, particle_n)
%F_ARCLENGTH

LS(1) = 0;
for i = 2:length(traj)
    for particle = 1:particle_n
        temp(particle) = traj(particle, i) - traj(particle, i - 1);
        
    end
    S(i) = sqrt(sum(temp.^2));
    LS(i) = LS(i-1) + S(i);
end

for i = 1:length(traj)
    VS(i) = 0;
    for particle = 1:particle_n
        VS(i) = VS(i) + 0.25 * (traj(particle, i)^2 + a)^2;
        for particleB = (particle + 1):particle_n
            VS(i) = VS(i) + eta/abs(traj(particle, i) - traj(particleB, i));
        end
    end
end
chiS = 0;
end

