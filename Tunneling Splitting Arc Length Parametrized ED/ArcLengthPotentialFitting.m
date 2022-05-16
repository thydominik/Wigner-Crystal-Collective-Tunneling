clc
clear all

% Constructing the arc length parametrized V(S) potential, then apply some
% fitting to that potential.

% CONSTANTS:---------------------------------------------------------------
eta = 20;   %dimensionless coulomb interaction strength

% -------------------------------------------------------------------------

% Addpath the files:
addpath("D:\BME PhD\Wigner Crystal Collective Tunneling\Data\Trajectories\eta 20\3 Particles");

%Script:
    % 1 - Load trajectories
    % 2 - Calculate the arc length parametrization {S(t); V(S)}
    % 3 - Apply various fitting to V(S) [quartic, quadratic, odd powers?]
FittingParameters = [];
q = 1;
for ind = 1:5:71

    %creating the name for each trajectories that will be loaded in
    name    = ind;
    nameSTR = ['P_' num2str(name)];

    %loading a trajectory that is previoulsy determined by a MC simulatiton
    trajectory_load = load(nameSTR);
    trajectory      = trajectory_load.position;

    %The equilibrium positons can be loaded in, or taken from the
    %trajectories end points, by constrction they will be the same.
    equilibrium_positions   = load('EqPos_eta20_alpha_5_20.mat');    % 4th column is the alpha value!    I know that this is called in every loop, it's fine for now.
    state                   = ind;                               % just to keep track of the indexing meaning: state = alpha
    eq_pos                  = equilibrium_positions.eqpos(:, state);  % eq_pos holds the equilibrium poisitions for a psecific alpha only
    alpha                   = eq_pos(end); disp(['Alpha= ', num2str(alpha)])

    %it's useful not to set in stone the division and particle number for
    %later cases
    [particle_n, N_division]    = size(trajectory);

    % The trajectories are calculated using this z time:
    z_time          = linspace(-1, 1, N_division);

    % Using the trajectories, creating the arc length param and the arc
    % length parametrized V(S) potential.

    [chiS, S, VS]   = f_arclength(trajectory, alpha, eta, z_time);
    S               = S - max(S)/2;
    VS              = VS - min(VS);

    % Fitting Method 1: a/2 * x^2 + b/4 * x^4 + c
    [gof1, fc1] = f_fitting_VS_1(S, VS);
    
    % Fitting Method 2: 0.25b * (x^2 - min(S)^2)^2
    [gof2, fc2] = f_fitting_VS_2(S, VS);

    % Fitting Method 3: 0.25b * (x - x0)^2 * (x + x0)^2 + c
    [gof3, fc3] = f_fitting_VS_3(S, VS);


    figure(1)
    clf(figure(1))
    hold on
    title(['S(z) at \alpha = ' num2str(alpha)])
    xlabel('z')
    ylabel('S(z)')
    plot(z_time, S, '.-')
    hold off

    figure(2)
    clf(figure(2))
    hold on
    title(['V(S) at \alpha = ' num2str(alpha)])
    xlabel('S')
    ylabel('V(S)')
    plot(S, VS, '.')
    hold off


    figure(3)
    clf(figure(3))
    hold on
    x = linspace(-1.3 * max(S), max(S) * 1.3, 300);
    VSf1 = 0.5 * fc1.a * x.^2 + 0.25 * fc1.b * x.^4 + fc1.c;
    VSf2 = 0.25 * fc2.b * (x.^2 - min(S)^2).^2 ; 
    VSf3 = 0.25 * fc3.b * (x - min(S)).^2 .* (x + min(S)).^2 + max(VS) - (0.25 * fc3.b * min(S)^4);
    disp(mean(abs(VSf1 - VSf2)))
    disp(mean(abs(VSf1 - VSf3)))
    plot(S, VS, 'o', 'DisplayName', 'Original V(S)')
    plot(x, VSf1,  'DisplayName', 'Method 1')
    plot(x, VSf2,'DisplayName', 'Method 2')
    plot(x, VSf3, 'DisplayName', 'Method 3')
    legend
    grid on
    hold off
    
    scaling(q, 1) = alpha;
    scaling(q, 2) = min(S)^2 * fc2.b;
    scaling(q, 3) = fc2.b;
    q = q + 1;
    fc2.b - ((1/min(S)^4) * 4 * max(VS))
    pause
end
%%
figure(4)
clf(figure(4))
hold on
a = scaling(:, 1);
b = scaling(:, 2);
c = scaling(:, 3);
plot(nonzeros(scaling(: ,1)), -nonzeros(scaling(:, 2)))
plot(nonzeros(scaling(: ,1)), nonzeros(scaling(:, 3)))
grid on
hold off