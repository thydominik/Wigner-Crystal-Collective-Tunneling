function [traj] = OneParticleStep(traj, ParticleNum, sigma, Nt)
%

RndInd = randi( Nt - 2) + 1;

Deviation                   = normrnd(traj(ParticleNum, RndInd), sigma);
traj(ParticleNum, RndInd)   = Deviation;

end

