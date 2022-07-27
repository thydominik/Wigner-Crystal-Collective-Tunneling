clc
clear all

%1 particle trajectory calculation with simulated annealing

%constants:

NoP         = 200;                  % # of points in a curve
AlphaJumps  = 0.5;                 % Increments in the alpha
AlphaValues = 1:AlphaJumps:10;      % Potential barrier parameter

for k = 1:length(AlphaValues)
    R(k) = 3/sqrt(AlphaValues(k));
end
%the time rescaling parameter set to some arbitrary number O(1)

%potential paramters:
E_0 = 0.71738;          % the energy scale of the hamiltonian
E_p = 0.71738;          % the energy scale of the potential
l_d = 161.07;           % Length scale of the potential

% Imaginary Time: 
z   = linspace(-1, 1, NoP);
dz  = abs(z(1) - z(2));

for stateInd = 1:length(AlphaValues)
    tic
    
    r       = R(stateInd);
    alpha   = AlphaValues(stateInd)

    iter    = 2*10^7;                 % # of MC iterations 
    EqPosIn = -sqrt(alpha);
    EqPosOut = sqrt(alpha);

    % Simulated temperature & Sigma (variance)
    T_init  = 10;
    T       = T_init * exp(-(linspace(0,40,iter))/2);
    sigma   = 0.1 * sqrt(T);

    [Position, Shift] = InitPos1Particle(alpha, NoP);
    CurrentAction = ActionCalc(Position, r, NoP, z, dz, alpha);

    PosAnalytic = sqrt(alpha)*tanh(atanh(z)*r*sqrt(alpha)/sqrt(2));
    ActionAnalytic = ActionCalc(PosAnalytic, r, NoP, z, dz, alpha);

    for i = 1:iter
        pos_new = Newstep(Position,NoP,sigma(i), alpha, z);
        E_new   = ActionCalc(pos_new, r, NoP, z, dz, alpha);
        E_diff  = CurrentAction - E_new;

        if E_diff > 0
            Position        = pos_new;
            CurrentAction   = E_new;
        elseif rand() <= exp(E_diff/T(i))
            Position = pos_new;
            CurrentAction = E_new;
        else
            discarded(i) = i;
        end

        E(i) = CurrentAction;
        if rem(i,100000) == 0
            %disp(i)
            disp(E(i) - ActionAnalytic)
            figure(2)
            clf(figure(2))
            hold on
            ylim([(-sqrt(alpha)-0.1) (sqrt(alpha)+0.1)])
            scatter(z,Position)
            plot(z, sqrt(alpha)*tanh(atanh(z)*r*sqrt(alpha)/sqrt(2)))
            xlabel(['i z(r); ', 'r = ',num2str(r),'; \alpha = ',num2str(alpha),'; E_p = ',num2str(E_p)],'FontSize',19)
            ylabel('\chi(z)','FontSize',19)
            hold off
            disp(['iter: ', num2str(i)])
        end

    end

    toc
    NameString = ['Traj_1p_' num2str(state)];

    IterData.NameString                 = NameString;
    IterData.NumberOfParticles          = 1;
    IterData.AlphaIncrements            = AlphaJumps;
    IterData.AlphaValues                = AlphaValues;
    IterData.NumberOfPoinsInTrajectory  = NoP;
    IterData.time                       = z;
    IterData.dt                         = dz;
    IterData.EquilibriumPositions       = [EqPosIn EqPosOut];
    IterData.RValue                     = r;
    IterData.IterationNumber            = iter;
    IterData.EnergySeries               = Energy(1:10:end);
    IterData.EnergyShift                = Shift;
    IterData.Trajectories               = Position;
    IterData.Action                     = E(end);
    IterData.AnalyticalAction           = ActionAnalytic;

    save(NameString, "IterData");
end
