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

