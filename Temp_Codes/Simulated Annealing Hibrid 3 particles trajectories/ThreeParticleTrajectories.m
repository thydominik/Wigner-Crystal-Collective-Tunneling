clc
clear all

FS = 1;     %figure switch 0 - figures off; 1 - figures on

% Initialize the problem
    % Get particle number, physical constants, R values,
NoP = 3;    %Number of Particles

%The alpha values:
alpha_jumps = 0.25;
alphavals       = 5 : alpha_jumps : 10;

eta = 20;                   % Dimensionless coulomb interaction strength
Nt  = 200;                  % Number of points in the trajectories
z   = linspace(-1, 1, Nt);  % Time variable
dz  = z(2) - z(1);          % Difference between two time steps


disp(['Number of Particles:                     ' num2str(NoP)])
disp(['dalpha:                                  ' num2str(alpha_jumps)])
disp(['\alpha_min: ' num2str(alphavals(1)) ', \alpha_max: ' num2str(alphavals(end))])
disp(['Number of points in one trajectory:      ' num2str(Nt)])
disp(['dz:                                      ' num2str(dz)])
disp(['Dimless Coulomb interaction parameter:   ' num2str(eta)])

% Equilibrium positions

EqPos = [];

for a = 1:length(alphavals)
    Potential = @(x) 0.25 * (x(1)^2 - alphavals(a))^2 + 0.25 * (x(2)^2 - alphavals(a))^2 + 0.25 * (x(3)^2 - alphavals(a))^2 + eta/abs(x(1) - x(2)) + eta/abs(x(1) - x(3)) + eta/abs(x(2) - x(3));

    options = optimset('TolFun', 1e-12, 'TolX', 1e-12, 'MaxFunEvals', 10^8, 'MaxIter', 10^8);
    x_start = [-sqrt(a)-1 -1 sqrt(a)];
    [x0, fval0] = fminsearch(Potential, x_start, options);
    EqPos(:, a) = x0;
    FuncVal = fval0;
end

if FS
    figure(1)
    clf(figure(1))
    hold on
    title('Equilibrium positions')
    xlabel('|\alpha|')
    ylabel('\chi_i')
    for i = 1:NoP
        plot(alphavals, EqPos(i, :), '.-', 'DisplayName', ['Particle ' num2str(i)])
    end
    y = linspace(4.5, max(alphavals));
    plot(y, -sqrt(y), 'k', 'DisplayName', '-sqrt(\alpha)')
    plot(y, sqrt(y), 'k', 'DisplayName', 'sqrt(\alpha)')
    legend
    grid on
    hold off
end

disp('Equilibrium positions calculated. --------------------------------------')

% Caclculate Trajectories -- Simulated Annealing

%R = linspace(2, 0.3, length(alphavals));   % 1 - 17
R = linspace(4, 0.8, length(alphavals));      % 18 - end
Exclude = 0.2;
for state = 1:length(alphavals)
    tic
    r       = R(state);   % Time scaling parameter
    alpha   = alphavals(state); %Specific state's alpha value

    iter = 4 * 10^6;    %# of iterations for one specific alpha

    for q = 1:NoP
        p_in(q)     = EqPos(q, state); 
        p_out(q)    = -EqPos(NoP - q + 1, state);
    end

    % Initializing the starting positions & determining the energy shift:
    [Position, Shift]   = InitPos3Particle(NoP, Nt, p_in, p_out, eta, alpha);
    CurrentAction       = ActionCalc(Position, r, alpha, eta, NoP, Nt, z, dz, Shift);

    % Simulated temperature & Sigma (variance)
    T_init = 20;
    T = T_init * exp(-(linspace(0,30,iter)));
    %T = [T linspace(T(end), 10^-16, iter/2)];

    sigma = 0.5 * sqrt(T);

    % Keeping track of the discarded moves & Energy/iteration:
    discarded = zeros(iter,1);
    Energy = zeros(iter,1);

    disp(["\alpha = " num2str(alpha)])
    disp("E_0 = " + num2str(CurrentAction, 15))
    disp("Shift = " + num2str(Shift))

    % SIMULATED ANNEALING PART:
    for i = 1:iter
        pos_new     = HibridStep3Part(Position, alpha, eta, NoP, Nt, sigma(i), Exclude, p_in, p_out, z);
        E_new       = ActionCalc(pos_new, r, alpha, eta, NoP, Nt, z, dz, Shift);
        E_diff      = CurrentAction - E_new;

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

        if rem(i,20000) == 0
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
            xline([1-Exclude (-1 + Exclude)])
            grid on
            hold off
            disp("iter= " + num2str(i) + "   "+ "E_0= " + num2str(CurrentAction, 10) + ' -- ' + num2str(round(i/iter * 100, 3)) + '%')
        end 

    end
    Action(2, state) = Energy(end);
    NameString = ['New_Traj_3p_' num2str(state)];
    save(NameString, 'Position')
    toc

    tic
    for ParticleIdx = 1:3
        iter2 = 10^6;
        T_init = 10^-5;
        T = T_init * exp(-(linspace(0,30,iter2)));
        sigma2 = 0.05 * sqrt(T);

        for IterIdx = 1:iter2
            pos_new     = OneParticleStep(Position, ParticleIdx, sigma2(IterIdx), Nt);
            E_new       = ActionCalc(pos_new, r, alpha, eta, NoP, Nt, z, dz, Shift);
            E_diff      = CurrentAction - E_new;

            if E_diff > 0
                Position        = pos_new;
                CurrentAction   = E_new;
            elseif rand() <= exp(E_diff/(0.2 * T(IterIdx)))
                Position        = pos_new;
                CurrentAction   = E_new;
            else
                discarded(end + 1) = IterIdx;
            end

            if rem(IterIdx, 20000) == 0
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
                xline([1-Exclude (-1 + Exclude)])
                grid on
                hold off
                disp("iter= " + num2str(IterIdx) + "   "+ "E_0= " + num2str(CurrentAction, 10))
            end
        end

    end
    toc
end





