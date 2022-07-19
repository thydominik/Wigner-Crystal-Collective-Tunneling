function [traj] = OneParticleStep(traj, ParticleNum, sigma, Nt)
    %
    
    RndInd = randi(Nt/2 - 1) + 1;
    
    Deviation                   = randn(1) * sigma + traj(ParticleNum, RndInd);
    
    traj(ParticleNum, RndInd)           = Deviation;
    traj(5 - ParticleNum + 1, Nt - RndInd + 1)   = -Deviation;
    
    while traj(ParticleNum, RndInd) < traj(ParticleNum, 1) || traj(ParticleNum, RndInd) > traj(ParticleNum, Nt/2 + 1)
        Deviation                                   = randn(1) * sigma + traj(ParticleNum, RndInd);
        traj(ParticleNum, RndInd)                   = Deviation;
        traj(5 - ParticleNum + 1, Nt - RndInd + 1)  = -Deviation;
    end

end

