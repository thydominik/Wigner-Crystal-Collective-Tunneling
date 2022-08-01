function [pos] = HibridStep5Part(pos, alpha, eta, NoP, Nt, sigma, Exclude, p_in, p_out, z)
    %HIBRIDSTEP5PART: Makes a hibrid step from pos, and gives back pos_new as a
    %new set of trajectories. pos - [array, double] set of trajs.; r - [double]
    %time scaling variable; alpha - [double] alpha parameter; eta - [double]
    %dimless Coulomb int. param.; NoP - [int] # of particles; Nt - [int] # of
    %points in trajectory; z - [array, double] time variable; dz - [double]
    %time step; sigma - [double] new step deviation from initial values; Exclude -
    % [double] excludes the derivative of the trajectories on z * exclude part.
    
    %Middle particle index
    Midx = floor(NoP/2) + 1;

    %Random points in the trajectories
    MiddleRndInd = randi(Nt/2 - 1) + 1;
    UpSideRndInd = randi(Nt/2 - (Nt/2 * Exclude)) + (Nt/2 * Exclude);
    DownSideRndInd = randi(Nt/2 - (Nt/2 * Exclude)) + (Nt/2 * Exclude);
    
    % First set the middle particle:
    
    DeviationMiddle                     = randn(1) * sigma + pos(Midx, MiddleRndInd); %normrnd(pos(Midx, MiddleRndInd), sigma);      %
    OriginalMiddleValue                 = pos(Midx, MiddleRndInd);
    pos(Midx, MiddleRndInd)             = DeviationMiddle;
    pos(Midx, Nt - MiddleRndInd + 1)    = -DeviationMiddle;
    
    while pos(Midx, MiddleRndInd) < p_in(Midx) || pos(Midx, MiddleRndInd) > 0 || pos(Midx, MiddleRndInd) > (z(MiddleRndInd) * abs(p_in(Midx)))
        DeviationMiddle = randn(1) * sigma + OriginalMiddleValue; %normrnd(OriginalMiddleValue, sigma);
        pos(Midx, MiddleRndInd)            = DeviationMiddle;
        pos(Midx, Nt - MiddleRndInd + 1)   = -DeviationMiddle;
    end  
    
    % Check if MiddleRndInd is in the exclusion range, if so calculate
    % equilibrium positions:
    if MiddleRndInd <= Nt/2 * Exclude
        Pm = pos(Midx, MiddleRndInd);      % Middle particle position, just a shorter notation
        Potential = @(x) 0.25 * (x(1)^2 - alpha)^2 + 0.25 * (x(2)^2 - alpha)^2 + 0.25 * (Pm^2 - alpha)^2 + 0.25 * (x(3)^2 - alpha)^2 + 0.25 * (x(4)^2 - alpha)^2 + eta * (1/abs(x(1) - x(2)) + 1/abs(x(1) - x(3)) + 1/abs(x(1) - x(4)) + 1/abs(x(1) - Pm) + 1/abs(x(2) - x(3)) + 1/abs(x(2) - x(4)) + 1/abs(x(2) - Pm) + 1/abs(x(3) - x(4)) + 1/abs(x(3) - Pm) + 1/abs(x(4) - Pm));
        
        options = optimset('TolFun', 1e-12, 'TolX', 1e-12, 'MaxFunEvals', 10^6, 'MaxIter', 10^6);
        x_start = [(-sqrt(alpha)-1) (-sqrt(alpha)) (sqrt(alpha)) (sqrt(alpha)+1)];
        [x0, fval0] = fminsearch(Potential, x_start, options);
        x0 = sort(x0);

        pos(1, MiddleRndInd) = x0(1);
        pos(2, MiddleRndInd) = x0(2);
        pos(4, MiddleRndInd) = x0(3);
        pos(5, MiddleRndInd) = x0(4);
        pos(5, Nt - MiddleRndInd + 1) = -x0(1);
        pos(4, Nt - MiddleRndInd + 1) = -x0(2);
        pos(2, Nt - MiddleRndInd + 1) = -x0(3);
        pos(1, Nt - MiddleRndInd + 1) = -x0(4);

    end

    %Now the side particles, the usual way

    DeviationSide1      = randn(1) * (0.8 * sigma) + pos(Midx - 1, DownSideRndInd); %normrnd(pos(Midx - 1, DownSideRndInd), sigma * 0.8);
    OriginalSide1Value  = pos(Midx - 1, DownSideRndInd);
    DeviationSide3      = randn(1) * (0.8 * sigma) + pos(Midx + 1, UpSideRndInd); %normrnd(pos(Midx + 1, UpSideRndInd), sigma * 0.8);
    OriginalSide3Value  = pos(Midx + 1, UpSideRndInd);

    pos(Midx - 1, DownSideRndInd)              = DeviationSide1;
    pos(Midx + 1, Nt - DownSideRndInd + 1)     = - DeviationSide1;

    while pos(Midx - 1, DownSideRndInd) < p_in(Midx - 1) || pos(Midx - 1, DownSideRndInd) > p_out(Midx - 1)
        D = randn(1) * sigma + OriginalSide1Value; %normrnd(OriginalSide1Value, sigma);
        pos(Midx - 1, DownSideRndInd)              = D;
        pos(Midx + 1, Nt - DownSideRndInd + 1)   = -D;
    end

    pos(Midx + 1, UpSideRndInd)            = DeviationSide3;
    pos(Midx - 1, Nt - UpSideRndInd + 1)   = - DeviationSide3;

    while pos(Midx + 1, UpSideRndInd) < p_in(Midx + 1) || pos(Midx + 1, UpSideRndInd) > p_out(Midx + 1)
        D = randn(1) * sigma + OriginalSide3Value; %normrnd(OriginalSide3Value, sigma);
        pos(Midx + 1, UpSideRndInd)             = D;
        pos(Midx - 1, Nt - UpSideRndInd + 1)    = -D;
    end
    
    Pm = pos(Midx, MiddleRndInd);      % Middle particle position, just a shorter notation
    Pd = pos(Midx - 1, DownSideRndInd);
    Pu = pos(Midx + 1, DownSideRndInd);

    Potential = @(x) 0.25 * (x(1)^2 - alpha)^2 + 0.25 * (x(2)^2 - alpha)^2 + 0.25 * (Pm^2 - alpha)^2 + 0.25 * (Pd^2 - alpha)^2 + 0.25 * (Pu^2 - alpha)^2 + eta * (1/abs(x(1) - Pd) + 1/abs(x(1) - Pm) + 1/abs(x(1) - Pu) + 1/abs(x(1) - x(2)) + 1/abs(Pd - Pm) + 1/abs(Pd - Pu) + 1/abs(Pd - x(2)) + 1/abs(Pm - Pu) + 1/abs(Pm - x(2)) + 1/abs(Pu - x(2)));

    options = optimset('TolFun', 1e-12, 'TolX', 1e-12, 'MaxFunEvals', 10^8, 'MaxIter', 10^8);
    x_start = [(-sqrt(alpha)-1) (sqrt(alpha)+1)];
    [x0, fval0] = fminsearch(Potential, x_start, options);
    x0 = sort(x0);

    pos(1, DownSideRndInd) = x0(1);
    pos(5, DownSideRndInd) = x0(2);
    pos(5, Nt - DownSideRndInd + 1) = -x0(1);
    pos(1, Nt - DownSideRndInd + 1) = -x0(2);

end

