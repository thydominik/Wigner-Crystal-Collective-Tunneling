%Simulated Annealing for general N particles. The calculation of the
%equilibrium positions are not included in this code, those needs to be
%calculated separately.

clc
clear all

%Loading the equlibrium positions:
eqpos   = load('EqPos_eta20_alpha_5_20_5p.mat');
eq_pos  = eqpos.eqpos;

%Getting the number of particles involved:
N = length(eq_pos(:, 1)) - 1;

%Trajectory selection:
TrajNum     = 5; %Parameter to skip a certain number of positions
TrajStart   = 50; %Starting from the TrajStart \alpha value
TrajFin     = length(eq_pos);

Action          = zeros(2, length(eq_pos)); %First column alpha; second column actions
Action(1, :)    = eq_pos(6, :);             %These are the alpha values

%Now this is sort of the tricky part: In order to get "good" trajectories,
%one should choose such r parameter, that ensures that the second point
%differs from hte trajecory's first point "significantly". This leads to
%less noise, and hopefully more accurete trajs. This needs to be saved
%alongside the trajectories!!!
R = linspace(1, 0.3, length(eq_pos)); %this is not yet tested
    %For now a contant is enough
R(:) = 1.5; 

%Now the Monte Carlo part:
for State = TrajStart:TrajNum:TrajFin
    tic
    %Parameters and constant of the simulation:
    r       = R(State);
    eps     = 10^-10;
    Nt      = 100;
    alpha   = eq_pos(N + 1, State);
    disp(["\alpha = " num2str(alpha)])
    
    rs      = 20;               %dimensionless interaction strength 18.813;
    iter    = 7 * 10^6;     %# of iterations

    %Time parameter:
    z_reduced   = linspace(-1 + eps, 1 - eps, Nt);
    z           = linspace(-1, 1, Nt);
    dz_reduced  = z_reduced(2) - z_reduced(1);
    dz          = z(2) - z(1);

    %Equilibrium positions:
    for q = 1:N
        p_in(q)     = eq_pos(q, State); 
        p_out(q)    = -eq_pos(N - q + 1, State);
    end

    %simulated temperature & Sigma (variance)
    T_init = 5;
    T = T_init * exp(-(linspace(0,30,iter)));

    sigma = 0.2 * T;

    %Initializing the starting positions & determining the energy shift:
    [Position, Shift]   = f_initpos(N, Nt, p_in, p_out, rs, alpha);
    E_0                 = f_actioncalc(Position, r, alpha, rs, N, Nt, z, dz, Shift);
    disp("E_0 = " + num2str(E_0, 15))
    disp("Shift = " + num2str(Shift))

    %Keeping track of the discarded moves & Energy/iteration:
    discarded = zeros(iter,1);
    E = zeros(iter,1);
    
    %Simulated Annealing part:
    for i = 1:iter
        pos_new     = f_newstep(Position, N, Nt, sigma(i), p_in, p_out);
        E_new       = f_actioncalc(pos_new, r, alpha, rs, N, Nt, z, dz, Shift);
        E_diff      = E_0 - E_new;

        if E_diff > 0
            Position    = pos_new;
            E_0         = E_new;
        elseif rand() <= exp(E_diff/T(i))
            Position    = pos_new;
            E_0         = E_new;
        else
            discarded(i) = i;
        end
        E(i) = E_0;

        if rem(i,5000) == 0
            figure(1)
            clf(figure(1))
            hold on
            plot(z, Position(1, :))
            plot(z, Position(2, :))
            plot(z, Position(3, :))
            plot(z, Position(4, :))
            plot(z, Position(5, :))
            plot(0.01 * 0.25 * (linspace(-5, 5, 100).^2 + alpha).^2 - 1, linspace(-5, 5, 100))
            xlim([-1 1])
            yline([p_in(1) p_in(2) p_in(3) p_in(4) p_in(5) 0])
            yline([p_out(1) p_out(2) p_out(3) p_out(4) p_out(5)])
            grid on
            hold off
            disp("iter= " + num2str(i) + "   "+ "E_0= " + num2str(E_0, 10))
        end   
    end
   
    Action(2, State) = E(end);
    NameString = ['Traj_5p_' num2str(State)];
    save(NameString, 'position')

    toc
end

