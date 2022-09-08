clc
clear all

format long

%% Structural Initialization:
FGS = 1;        % Figure switch
FE  = 1;        % ErrorFinding Figure switch
MS  = 1;        % Method switch: MS = {1, 2, 3} = {Standard MC, Restricted MC, No derivatives}
    % For Option two how many spwees do you want?
    if MS == 2
        SweepCount = 2;
    end

PN  = 7;        % Number of Particles

%% Constants of the calculation:
NoP             = 160;              % Number of points in each trajectory
AlphaJumps      = 0.5;              % Increments in Alpha
AlphaValues     = 11:AlphaJumps:14;  % 'Potential barrier' values
iter            = 1.4 * 10^8; %11.2 * 10^7;             % Number of MC iterations
R0 = 0.35 * 13;         % 0.25 is low for a = 11; but after 15 it could be good again.
for k = 1:length(AlphaValues)
    R(k) = R0 / sqrt(AlphaValues(k)); %the time rescaling parameter set to some arbitrary number O(1)
end

Eta = 20;  % Potential fitted value of the dimensionless Coulomb interaction 18.813...
z   = linspace(-1, 1, NoP);     % Imaginary time variable
dz  = z(2) - z(1);              % dt in imaginary time

disp(['Number of Particles:                     ' num2str(PN)])
disp(['dalpha:                                  ' num2str(AlphaJumps)])
disp(['\alpha_min: ' num2str(AlphaValues(1)) ', \alpha_max: ' num2str(AlphaValues(end))])
disp(['Number of points in one trajectory:      ' num2str(NoP)])
disp(['dz:                                      ' num2str(dz)])
disp(['Dimless Coulomb interaction parameter:   ' num2str(Eta)])

%% Calculating the Equilibrium positions

EqPos = [];

alpha = AlphaValues;
N = length(alpha);

EqPos = [];

for i = 1:N
    a = alpha(i);
    Potential = @(x) 0.25 * ((x(1)^2 - a)^2 + (x(2)^2 - a)^2 + (x(3)^2 - a)^2 + (x(4)^2 - a)^2 + (x(5)^2 - a)^2 + (x(6)^2 - a)^2 + (x(7)^2 - a)^2) + Eta * (1/abs(x(1) - x(2)) + 1/abs(x(1) - x(3)) + 1/abs(x(1) - x(4)) + 1/abs(x(1) - x(5)) + 1/abs(x(1) - x(6)) + 1/abs(x(1) - x(7)) + 1/abs(x(2) - x(3)) + 1/abs(x(2) - x(4)) + 1/abs(x(2) - x(5)) + 1/abs(x(2) - x(6)) + 1/abs(x(2) - x(7)) + 1/abs(x(3) - x(4)) + 1/abs(x(3) - x(5)) + 1/abs(x(3) - x(6)) + 1/abs(x(3) - x(7)) + 1/abs(x(4) - x(5)) + 1/abs(x(4) - x(6)) + 1/abs(x(4) - x(7)) + 1/abs(x(5) - x(6)) + 1/abs(x(5) - x(7)) + 1/abs(x(6) - x(7)));
    % options = optimset('Display','iter','PlotFcns',@optimplotfval, 'TolFun', 1e-8, 'TolX', 1e-8);
    options = optimset('TolFun', 1e-300, 'TolX', 1e-300, 'MaxFunEvals', 10^16, 'MaxIter', 10^16);
    x_start = [(-sqrt(a)-2) (-sqrt(a)) (-sqrt(a)+1) (-sqrt(a)+2.1) (sqrt(a)-2) (sqrt(a)-1) (sqrt(a)+1)];
    %x_start = [(-sqrt(a)-2) (-sqrt(a)) (-sqrt(a)+1) 0 (sqrt(a)-2) (sqrt(a)-1) (sqrt(a)+1)];
    [x0, fval0] = fminsearch(Potential, x_start, options);
    EqPos(i, :) = sort(x0);
    FuncVal = fval0;
end

disp('Done with 7 particle')

figure(1)
clf(figure(1))
hold on
title('Equilibrium positions for 7 particle')
xlabel('\alpha')
ylabel('\chi')
plot(alpha, EqPos(:, 1), 'o-')
plot(alpha, EqPos(:, 2), 'o-')
plot(alpha, EqPos(:, 3), 'o-')
plot(alpha, EqPos(:, 4), 'o-')
plot(alpha, EqPos(:, 5), 'o-')
plot(alpha, EqPos(:, 6), 'o-')
plot(alpha, EqPos(:, 7), 'o-')
plot(alpha, sqrt(alpha), 'k')
plot(alpha, -sqrt(alpha), 'k')
yline(0)

if FGS
    figure(1)
    clf(figure(1))
    hold on
    title('Equilibrium positions')
    xlabel('|\alpha|')
    ylabel('\chi_i')
    for i = 1:PN
        plot(AlphaValues, EqPos(:, i), '.-', 'DisplayName', ['Particle ' num2str(i)])
    end
    y = linspace(13, max(AlphaValues));
    plot(y, -sqrt(y), 'k', 'DisplayName', '-sqrt(\alpha)')
    plot(y, sqrt(y), 'k', 'DisplayName', 'sqrt(\alpha)')
    legend
    grid on
    hold off
end

disp('Equilibrium positions calculated. --------------------------------------')

%% Caclculate Trajectories -- Simulated Annealing
Exclude = 0;

for stateInd = 1:length(AlphaValues)
    % Simulated temperature & Sigma (variance)
    T_init = 40;
    T = T_init * exp(-(linspace(0,36,iter)));

    sigma = 0.01 * sqrt(T);

    Exclude = 0.0;

    r       = R(stateInd);   % Time scaling parameter
    alpha   = AlphaValues(stateInd); %Specific state's alpha value

    for q = 1:PN
        p_in(q)     = EqPos(stateInd, q);
        p_out(q)    = -EqPos(stateInd, PN - q + 1);
    end

    % Initializing the starting positions & determining the energy shift:
    [Position, Shift]   = InitPos7Particle(PN, NoP, p_in, p_out, Eta, alpha);

    if FGS == 1 && FE == 1
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

    CurrentAction   = ActionCalc(Position, r, alpha, Eta, PN, NoP, z, dz, Shift);

    % Keeping track of the discarded moves & Energy/iteration:
    discarded   = zeros(iter, 1);
    Energy      = zeros(iter, 1);

    disp(["\alpha = " num2str(alpha)])
    disp("E_0 = " + num2str(CurrentAction, 15))
    disp("Shift = " + num2str(Shift))

    %tic
    for i = 1:iter
        if MS == 1
            % Option 1: Basic 5 particle MC/ SA
            pos_new = NewStepStandard(Position, alpha, Eta, PN, NoP, sigma(i), p_in, p_out, z);
        elseif MS == 2
            % Option 2: Restricted Range 5 particle MC/SA
            pos_new = HibridStep5Part(Position, alpha, Eta, PN, NoP, sigma(i), Exclude, p_in, p_out, z);
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

        if FGS
            if rem(i, 20000) == 0
                figure(2)
                clf(figure(2))
                hold on
                x = linspace(-1, 1, NoP/2);
                f = linspace(p_in(4), p_out(4), NoP/2);
                plot(x, f, 'k')  
                plot(z, Position(1, :))
                plot(z, Position(2, :))
                plot(z, Position(3, :))
                plot(z, Position(4, :))
                plot(z, Position(5, :))
                plot(z, Position(6, :))
                plot(z, Position(7, :))
                plot(0.01 * 0.25 * (linspace(-6, 6, 100).^2 - alpha).^2 - 1, linspace(-5, 5, 100))
                xlim([-1 1])
                yline(p_in)
                yline(p_out)
                xline([1-Exclude (-1 + Exclude)])
                grid on
                hold off
                disp("iter= " + num2str(i) + "   "+ "E_0= " + num2str(CurrentAction, 20) + ' ' + num2str(i/iter * 100) + '%')
            end
        else
            if rem(i, 20000) == 0
                disp("iter= " + num2str(i) + "   "+ "E_0= " + num2str(CurrentAction, 10))
            end
        end

    end
    disp('Done with Primary MC')
   
    if MS == 2
        for sweepInd = 1:SweepCount
            for ParticleInd = 1:5
                iter2 = 5 * 10^6;

                T_init2 = 10^-5;
                T2 = T_init2 * exp(-(linspace(0,35,iter2)));
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

                    if rem(IterInd2, 20000) == 0 && FGS == 1
                        figure(4)
                        clf(figure(4))
                        hold on
                        plot(z, Position(1, :))
                        plot(z, Position(2, :))
                        plot(z, Position(3, :))
                        plot(z, Position(4, :))
                        plot(z, Position(5, :))
                        plot(0.01 * 0.25 * (linspace(-6, 6, 100).^2 - alpha).^2 - 1, linspace(-5, 5, 100))
                        xlim([-1 1])
                        yline([p_in(1) p_in(2) p_in(3) p_in(4) p_in(5)])
                        yline([p_out(1) p_out(2) p_out(3) p_out(4) p_out(5)])
                        xline([1-Exclude (-1 + Exclude)])
                        grid on
                        hold off
                        disp("iter= " + num2str(IterInd2) + "   "+ "E_0= " + num2str(CurrentAction, 10))
                    end
                end
            end
        end
    end

   % time(stateInd) = toc;

    NameString = ['Traj_7p_' num2str(stateInd)];
    save(NameString, 'Position')

    IterData.NameString                 = NameString;
    IterData.NumberOfParticles          = PN;
    IterData.AlphaIncrements            = AlphaJumps;
    IterData.AlphaValues                = AlphaValues;
    IterData.AlphaValue                 = alpha;
    IterData.CoulombStrength            = Eta;
    IterData.NumberOfPoinsInTrajectory  = NoP;
    IterData.time                       = z;
    IterData.dt                         = dz;
    IterData.EquilibriumPositions       = EqPos;
    IterData.RValues                    = r;
    %IterData.Excludedregion             = Exclude;
    IterData.IterationNumber            = iter;
    %IterData.TemperatureValues          = T;
    %IterData.Sigma                      = sigma;
    %IterData.EnergySeries               = Energy(1:10:end);
    %IterData.EnergyShift                = Shift;
    %IterData.OneParticleIteration       = iter2;
    %IterData.OneParticleTemperature     = T1;
    %IterData.AcceptancefactorForTemp    = 0.2;
    %IterData.OneParticleSigma           = sigma2;
    %IterData.DiscardedSteps             = discarded;
    %IterData.oneParticelEnergies        = E1(1:10:end);
    %IterData.ElapsedTime                = time(end);
    IterData.Trajectories               = Position;
    IterData.Action                     = Energy(end);
    IterData.MinAction                  = min(Energy);

    save(['TrajectoryData_7p_' num2str(stateInd)], "IterData");
end

