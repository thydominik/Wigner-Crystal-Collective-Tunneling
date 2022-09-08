clc
clear all
format long

%% Structural Initialization:
FS  = 1;        % Figure switch
FE  = 1;        % ErrorFinding Figure switch
MS  = 1;        % Method switch: MS = {1, 2, 3} = {Standard MC, Restricted MC, No derivatives}
    % For Option two how many spwees do you want?
    if MS == 2
        SweepCount = 2;
    end

PN  = 3;        % Number of Particles

%% Constants of the calculation:
NoP             = 200;              % Number of points in each trajectory
AlphaJumps      = 0.5;             % Increments in Alpha
AlphaValues     = 5.5:AlphaJumps:13;  % 'Potential barrier' values
iter            = 5 * 10^6;             % Number of MC iterations
R0 = 4;
for k = 1:length(AlphaValues)
    R(k) = R0 / sqrt(AlphaValues(k)); %the time rescaling parameter set to some arbitrary number O(1)
end

Eta = 20;  % Potential fitted value of the dimensionless Coulomb interaction 18.813...
z   = linspace(-1, 1, NoP);     % Imaginary time variable
dz  = z(2) - z(1);              % dt in imaginary time

%% Calculating the Equilibrium Positions:

EqPos = [];
warning('Always check the equilibrium positions when changing the AlphaJumps or the initial and final value of alpha in AlphaValues before running the whole calculation!')
for a = 1:length(AlphaValues)
    Potential = @(x) 0.25 * (x(1)^2 - AlphaValues(a))^2 + 0.25 * (x(2)^2 - AlphaValues(a))^2 + 0.25 * (x(3)^2 - AlphaValues(a))^2 + Eta/abs(x(1) - x(2)) + Eta/abs(x(1) - x(3)) + Eta/abs(x(2) - x(3));

    options = optimset('TolFun', 1e-14, 'TolX', 1e-14, 'MaxFunEvals', 10^9, 'MaxIter', 10^9);
    x_start = [-sqrt(a)-2 -1 sqrt(a)];
    [x0, fval0] = fminsearch(Potential, x_start, options);
    EqPos(:, a) = x0;
    FuncVal(a)  = fval0;
end

if FS
    figure(1)
    clf(figure(1))
    hold on
    title('Equilibrium positions')
    xlabel('|\alpha|')
    ylabel('\chi_i')
    for i = 1:PN
        plot(AlphaValues, EqPos(i, :), '.-', 'DisplayName', ['Particle ' num2str(i)])
    end
    y = linspace(4.5, max(AlphaValues));
    plot(y, -sqrt(y), 'k', 'DisplayName', '-sqrt(\alpha)')
    plot(y, sqrt(y), 'k', 'DisplayName', 'sqrt(\alpha)')
    legend
    grid on
    hold off
end

disp('Equilibrium positions calculated. --------------------------------------')

disp(['Number of Particles:                     ' num2str(PN)])
disp(['dalpha:                                  ' num2str(AlphaJumps)])
disp(['\alpha_min: ' num2str(AlphaValues(1)) ', \alpha_max: ' num2str(AlphaValues(end))])
disp(['Number of points in one trajectory:      ' num2str(NoP)])
disp(['dz:                                      ' num2str(dz)])
disp(['Dimless Coulomb interaction parameter:   ' num2str(Eta)])
%% Simulated Annealing:

for stateInd = 1:length(AlphaValues)
    % Simulated temperature & Sigma (variance)
    T_init  = 10;
    T       = T_init * exp(-(linspace(0, 34, iter))/1);
    sigma   = 0.1 * sqrt(T);

    Exclude = 0.2;

    r       = R(stateInd);              % smae as R but we need only 1 for a specific calculation
    alpha   = AlphaValues(stateInd);
    
    for q = 1:PN
        p_in(q)     = EqPos(q, stateInd); 
        p_out(q)    = -EqPos(PN - q + 1, stateInd);
    end

    % Initializing the starting positions:
    [Position, Shift]   = InitPos3Particle(PN, NoP, p_in, p_out, Eta, alpha);
    
    if FS == 1 && FE == 1
        figure(3)
        clf(figure(3))
        hold on
        title('Initial trajectory')
        xlabel('z')
        ylabel('\chi_i')
        for i = 1:PN
            plot(z, Position(i, :), '.-', 'DisplayName', ['Particle ' num2str(i)])
        end
        yline(0)
        xline(0)
        legend
        grid on
        hold off
    end

    CurrentAction = ActionCalc(Position, r, alpha, Eta, PN, NoP, z, dz, Shift);
    
    % Keeping track of the discarded moves & Energy/iteration:
    discarded   = zeros(iter,1);
    Energy      = zeros(iter,1);

    disp(["\alpha = " num2str(alpha)])
    disp("E_0 = " + num2str(CurrentAction, 15))
    disp("Shift = " + num2str(Shift))

    tic
    for i = 1:iter
        if MS == 1
            % Option 1: Basic 3 particle MC
            pos_new = NewStepStandard(Position, NoP, sigma(i), p_in(1), p_out(1), p_in(2), p_out(2), p_in(3), p_out(3), z, dz);
        elseif MS == 2   
            % Option 2: Restricted range 3 particle MC
            pos_new = NewStepRestricted(Position, alpha, Eta, NoP, sigma(i), Exclude, p_in, p_out, z);
        elseif MS == 3
            % Option 3: No derivatives on the sides MC
            pos_new = HibridStep5PartEasy(Position, alpha, Eta, PN, NoP, sigma(i), Exclude, p_in, p_out, z);
        end
        E_new   = ActionCalc(pos_new, r, alpha, Eta, PN, NoP, z, dz, Shift);
        E_diff  = CurrentAction - E_new;

        if E_diff > 0
            Position        = pos_new;
            CurrentAction   = E_new;
        elseif rand() <= exp(E_diff/T(i))
            Position        = pos_new;
            CurrentAction   = E_new;
        else
            discarded(i) = i;
        end

        Energy(i) = CurrentAction;

        if rem(i, 20000) == 0 && FS == 1
            figure(2)
            clf(figure(2))
            hold on
            plot(z, Position(1, :))
            plot(z, Position(2, :))
            plot(z, Position(3, :))
            plot(0.01 * 0.25 * (linspace(-6, 6, 100).^2 - alpha).^2 - 1, linspace(-5, 5, 100))
            xlim([-1 1])
            yline([p_in(1) p_in(2) p_in(3)])
            yline([p_out(1) p_out(2) p_out(3)])
            plot(linspace(-1, 1, 25), linspace(Position(2, 1), Position(2, end), 25))
            grid on
            hold off
            disp("iter= " + num2str(i) + "   "+ "E_0= " + num2str(CurrentAction, 10) + ' -- ' + num2str(round(i/iter * 100, 3)) + '%')
        end
    end

    if MS == 2
        for sweepInd = 1:SweepCount
            for ParticleInd = 1:3
                iter2 = 10^6;

                T_init2 = 10^-5;
                T2 = T_init2 * exp(-(linspace(0,40,iter2)));
                sigma2 = 0.05 * sqrt(T2);

                for IterInd2 = 1:iter2
                    pos_new     = OneParticleStep(Position, ParticleInd, sigma2(IterInd2), NoP);
                    E_new       = ActionCalc(pos_new, r, alpha, Eta, PN, NoP, z, dz, Shift);
                    E_diff      = CurrentAction - E_new;

                    if E_diff > 0
                        Position        = pos_new;
                        CurrentAction   = E_new;
                    elseif rand() <= exp(E_diff/(0.2 * T2(IterInd2)))
                        Position        = pos_new;
                        CurrentAction   = E_new;
                    else
                        discarded(end + 1) = IterInd2;
                    end

                    if rem(IterInd2, 20000) == 0 && FS == 1
                        figure(4)
                        clf(figure(4))
                        hold on
                        plot(z, Position(1, :))
                        plot(z, Position(2, :))
                        plot(z, Position(3, :))
                        plot(0.01 * 0.25 * (linspace(-6, 6, 100).^2 - alpha).^2 - 1, linspace(-5, 5, 100))
                        xlim([-1 1])
                        yline([p_in(1) p_in(2) p_in(3)])
                        yline([p_out(1) p_out(2) p_out(3)])
                        xline([1-Exclude (-1 + Exclude)])
                        grid on
                        hold off
                        disp("iter= " + num2str(IterInd2) + "   "+ "E_0= " + num2str(CurrentAction, 10))
                    end
                end
            end
        end
    end

    %
    [gofmtx] = TanFitting(Position, NoP, PN, z);
    FittedPosition = zeros(3, NoP);
    for particleInd = 1:3
        if particleInd == 1 || particleInd == 3
            a = gofmtx(particleInd, 1);
            b = gofmtx(particleInd, 2);
            c = gofmtx(particleInd, 3);
            d = gofmtx(particleInd, 4);
            FittedPosition(particleInd, :) = a + b.*tanh(atanh(z).*c + d);
        else
            FittedPosition(particleInd, :) = abs(Position(2, 1)) .* tanh(atanh(z).*gofmtx(particleInd, 1));
        end
    end
    %
    if FS == 1
        figure(5)
        clf(figure(5))
        hold on
        plot(z, FittedPosition(1, :), 'o')
        plot(z, FittedPosition(2, :), 'o')
        plot(z, FittedPosition(3, :), 'o')
        plot(z, Position(1, :))
        plot(z, Position(2, :))
        plot(z, Position(3, :))
        grid on
        hold off
    end
    disp('fitted action: ')
    num2str(ActionCalc(FittedPosition, r, alpha, Eta, PN, NoP, z, dz, Shift))
    disp('MC action: ')
    num2str(Energy(end))

    time(stateInd) = toc;
    if MS == 1
        NameString = ['Traj_1p_STDMC' num2str(stateInd)];
    elseif MS == 2
        NameString = ['Traj_1p_RESTMC' num2str(stateInd)];
        IterData.SweepCount = SweepCount;
        IterData.Excluded   = Exclude;
    elseif MS == 3
        NameString = ['Traj_1p_SIMPLEMC' num2str(stateInd)];
    end
    

    IterData.NameString                 = NameString;
    IterData.NumberOfParticles          = PN;
    IterData.AlphaIncrements            = AlphaJumps;
    IterData.AlphaValues                = AlphaValues;
    IterData.NumberOfPoinsInTrajectory  = NoP;
    IterData.time                       = z;
    IterData.dt                         = dz;
    IterData.EquilibriumPositions       = EqPos;
    IterData.RValue                     = r;
    IterData.IterationNumber            = iter;
    %IterData.EnergySeries               = E(1:10:end);
    IterData.EnergyShift                = Shift;
    IterData.Trajectories               = Position;
    IterData.Action                     = Energy(end);
    IterData.AnalyticalAction           = 0;
    IterData.FittedAction               = ActionCalc(FittedPosition, r, alpha, Eta, PN, NoP, z, dz, Shift);
    IterData.FittedTrajectories         = FittedPosition;
    save(NameString, "IterData");

end



