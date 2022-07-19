function [pos] = HibridStep3Part(pos, alpha, eta, NoP, Nt, sigma, Exclude, p_in, p_out, z)
    %HIBRIDSTEP3PART: Makes a hibrid step from pos, and gives back pos_new as a
    %new set of trajectories. pos - [array, double] set of trajs.; r - [double]
    %time scaling variable; alpha - [double] alpha parameter; eta - [double]
    %dimless Coulomb int. param.; NoP - [int] # of particles; Nt - [int] # of
    %points in trajectory; z - [array, double] time variable; dz - [double]
    %time step; sigma - [double] new step deviation from initial values; Exclude -
    % [double] excludes the derivative of the trajectories on z * exclude part.
    
    MiddleRndInd = randi(Nt/2 - 1) + 1;
    UpSideRndInd = randi(Nt/2 - (Nt/2 * Exclude)) + (Nt/2 * Exclude);
    DownSideRndInd = randi(Nt/2 - (Nt/2 * Exclude)) + (Nt/2 * Exclude);
    
    % First set te middle particle:
    
    DeviationMiddle = normrnd(pos(2, MiddleRndInd), sigma);      %
    OriginalMiddleValue = pos(2, MiddleRndInd);
    pos(2, MiddleRndInd)            = DeviationMiddle;
    pos(2, Nt - MiddleRndInd + 1)   = -DeviationMiddle;
    
    while pos(2, MiddleRndInd) < p_in(2) || pos(2, MiddleRndInd) > 0 || pos(2, MiddleRndInd) > (z(MiddleRndInd) * abs(p_in(2)))
        DeviationMiddle = normrnd(OriginalMiddleValue, sigma);
        pos(2, MiddleRndInd)            = DeviationMiddle;
        pos(2, Nt - MiddleRndInd + 1)   = -DeviationMiddle;
    end  
    
    % Check if MiddleRndInd is in the exclusion range, if so calculate
    % equilibrium positions:
    if MiddleRndInd <= Nt/2 * Exclude
        Pm = pos(2, MiddleRndInd);      % Middle particle position, just a shorter notation
        Potential = @(x) 0.25 * (x(1)^2 - alpha)^2 + 0.25 * (Pm^2 - alpha)^2 + 0.25 * (x(2)^2 - alpha)^2 + eta/abs(x(1) - x(2)) + eta/abs(x(1) - Pm) + eta/abs(Pm - x(2));
        
        options = optimset('TolFun', 1e-12, 'TolX', 1e-12, 'MaxFunEvals', 10^6, 'MaxIter', 10^6);
        x_start = [-sqrt(alpha) sqrt(alpha)];
        [x0, fval0] = fminsearch(Potential, x_start, options);
        x0 = sort(x0);

        pos(1, MiddleRndInd) = x0(1);
        pos(3, Nt - MiddleRndInd + 1) = -x0(1);
        pos(3, MiddleRndInd) = x0(2);
        pos(1, Nt - MiddleRndInd + 1) = -x0(2);

    end

    %Now the side particles, the usual way

    DeviationSide1      = normrnd(pos(1, DownSideRndInd), sigma * 0.8);
    OriginalSide1Value  = pos(1, DownSideRndInd);
    DeviationSide3      = normrnd(pos(3, UpSideRndInd), sigma * 0.8);
    OriginalSide3Value  = pos(3, UpSideRndInd);

    pos(1, DownSideRndInd)              = DeviationSide1;
    pos(3, Nt - DownSideRndInd + 1)     = - DeviationSide1;

    while pos(1, DownSideRndInd) < p_in(1) || pos(1, DownSideRndInd) > p_out(1)
        D = normrnd(OriginalSide1Value, sigma);
        pos(1, DownSideRndInd)              = D;
        pos(3, Nt - DownSideRndInd + 1)   = -D;
    end

    pos(3, UpSideRndInd)            = DeviationSide3;
    pos(1, Nt - UpSideRndInd + 1)   = - DeviationSide3;

    while pos(3, UpSideRndInd) < p_in(3) || pos(3, UpSideRndInd) > p_out(3)
        D = normrnd(OriginalSide3Value, sigma);
        pos(3, UpSideRndInd)            = D;
        pos(1, Nt - UpSideRndInd + 1) = -D;
    end


end

