clc
clear all

%Calculating the tunneling splitting using an effective 1 particle arc
%length paramerized ED calculation:

% CONSTANTS:----------------------------------------------------------------------
eta = 20;   %dimensionless coulomb interaction strength

% --------------------------------------------------------------------------------
addpath("D:\BME PhD\Wigner Crystal Collective Tunneling\Data\Trajectories\eta 20\5 particles");
addpath("D:\BME PhD\Wigner Crystal Collective Tunneling\Data\Trajectories\eta 20\3 Particles");
%Script:
    % 1 - Load trajectories
    % 2 - Calculate the arc length parametrization {S(t); V(S)}
    % 3 - Apply various fitting to V(S) [quartic, quadratic, odd powers?]
    % 4 - Solve Schrödinger with ED
new_beta = [];
new_alpha = [];
alphaa = [];
for a_ind = 1:5:71      %Loop for a set of trajectories, this might have to be changes if theres is another set of trajectories

    %creating the name for each trajectories that will be loaded in
    name    = a_ind;
    % nameSTR = ['Traj_5p_' num2str(name)];
    nameSTR = ['P_' num2str(name)];

    %loading a trajectory that is previoulsy determined by a MC simulatiton
    trajectory_load = load(nameSTR);
    trajectory      = trajectory_load.position;

    %The equilibrium positons can be loaded in, or taken from the
    %trajectories end points, by constrction they will be the same.
    equilibrium_positions   = load('EqPos_eta20_alpha_5_20.mat');    % 4th column is the alpha value!    I know that this is called in every loop, it's fine for now.
    state                   = a_ind; %(a_ind - 1) /5 - 1;                                % just to keep track of the indexing meaning: state = alpha
    eq_pos                  = equilibrium_positions.eqpos(:, state);  % eq_pos holds the equilibrium poisitions for a psecific alpha only
    alpha                   = eq_pos(end); disp(['Alpha= ', num2str(alpha)])

    %it's useful not to set in stone the division and particle number for
    %later cases
    [particle_n, N_division]    = size(trajectory);
    
    % The trajectories are calculated using this z time:
    z_time          = linspace(-1, 1, N_division);

    % Using the trajectories, creating the arc length param and the arc
    % length parametrized V(S) potential.
    figure(6)
    clf(figure(6))
    hold on
    plot(z_time, trajectory(3, :))
    hold off

    [chiS, S, VS]   = f_arclength(trajectory, alpha, eta, z_time);
    S               = S - max(S)/2;
    VS              = VS - min(VS);
    
    % Now to have an even distribution of points, I do a slpine
    % interpolation on the V(S) data:
    N_interp    = 1000; %round((2 * max(S)) / (S(2) - S(1)))   % # of points in the interpolation
    S_interp    = linspace(S(1), S(end), N_interp); % Points ióon the arc length scale
    dS          = S_interp(2) - S_interp(1);
    VS_interp   = interp1(S, VS, S_interp, 'Spline');
    
    %eigs(zeros(round(N_interp * 1.5)), 10, 'smallestreal');

    % Fitting the V(S) potential with a quadratic potential
    reach = round(N_interp * (32/50));
    [gof2, fc2] = f_fitting_VS_2(S_interp(1:reach), VS_interp(1:reach));

    % Fitting the V(S) potential with a quartic potential
    [gof1, fc1] = f_fitting_VS_1(S_interp(1:reach), VS_interp(1:reach));

    new_beta(a_ind) = fc1.b;
    new_alpha(a_ind) = fc1.a;
    alphaa(a_ind) = alpha;
    % Now we have to build up the Potentials:
    S_min = (-1) * N_interp * 0.5 * dS + min(S_interp);
    k = 1;
    for i = 1:(2*N_interp)
        S_new = S_min + dS * (i - 1);
        % 1ST FOR THE QUADRATIC:
        if i <= N_interp/2
            Fitted_VS_2(i) = 0.25 * fc2.b * (S_new.^2 - min(S)^2).^2;
        elseif i > N_interp/2 && i <= (3/2 * N_interp)
            Fitted_VS_2(i) = VS_interp(k);
        else
            Fitted_VS_2(i) = 0.25 * fc2.b * (S_new.^2 - min(S)^2).^2;
        end
        Fitted_S_2(i) = S_new;

        % 2ND FOR THE QUARTIC:
        if i <= N_interp/2
            Fitted_VS_4(i) = 0.25 * fc1.b * S_new.^4 + 0.5 * fc1.a * S_new.^2 + fc1.c;
        elseif i > N_interp/2 && i <= (3/2 * N_interp)
            Fitted_VS_4(i) = VS_interp(k);
            k = k + 1;
        else
            Fitted_VS_4(i) = 0.25 * fc1.b * S_new.^4 + 0.5 * fc1.a * S_new.^2 + fc1.c;
        end
        Fitted_S_4(i) = S_new;
    end

    [Psi2, Spectra2] = Schrodinger_VSF(Fitted_S_2, Fitted_VS_2);
    dE2(a_ind) = Spectra2(2) - Spectra2(1);
    [Psi4, Spectra4] = Schrodinger_VSF(Fitted_S_4, Fitted_VS_4);
    dE4(a_ind) = Spectra4(2) - Spectra4(1);
    aaa(a_ind) = alpha;


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
    %plot(S_interp, VS_interp, 'o')
    hold off

    figure(3)
    clf(figure(3))
    hold on
    x = linspace(min(S) - abs(min(S)/2), min(S) + abs(min(S)/2), 100);
    title(['fit V(S) with 2nd ord. at \alpha = ' num2str(alpha)])
    xlabel('S')
    ylabel('V(S)')
    plot(x, fc2.b * (x - min(S)).^2, '.')
    plot(S_interp, VS_interp, 'o')
    plot(S_interp(1:reach), VS_interp(1:reach), 'o')
    hold off

    figure(4)
    clf(figure(4))
    hold on
    x = linspace(min(S) - abs(min(S)/2), min(S) + abs(min(S)/2), 100);
    title(['fit V(S) with 4th ord. at \alpha = ' num2str(alpha)])
    xlabel('S')
    ylabel('V(S)')
    plot(x, 0.5 * fc1.a * x.^2 + 0.25 * fc1.b * x.^4 + fc1.c, '.')
    plot(S_interp, VS_interp, 'o')
    plot(S_interp(1:reach), VS_interp(1:reach), 'o')
    hold off

    figure(5)
    clf(figure(5))
    hold on
    title(['fitted V(S) with 2nd ord. at \alpha = ' num2str(alpha)])
    xlabel('S')
    ylabel('V(S)')
    plot(Fitted_S_2, Fitted_VS_2, '.-', 'DisplayName', 'quadratic')
    plot(Fitted_S_4, Fitted_VS_4, '.-', 'DisplayName', 'quartic')
    plot(S_interp, VS_interp, 'o', 'DisplayName', 'Interpolated og')
    xline([min(S), max(S)], 'DisplayName', '\pm S_0')
    legend
    hold off

    

end

%%

M_data = load('E_Schrodinger_3e_eta_20.00_beta_0.01_N_100.dat');
data = load("data.mat");
data = data.data;

figure(6)
    clf(figure(6))
    hold on
    %title(['\Delta_2 and \Delta_4 at \alpha = ' num2str(alpha)])
    xlabel('|\alpha|', 'FontSize', 20)
    ylabel('\Delta', 'FontSize', 20)
    plot(abs(nonzeros(aaa)) - 4.45, nonzeros(dE2) .* data(:, 2), '.-', 'DisplayName', 'quadratic fit for eff. V(S) - ED')
    plot(abs(nonzeros(aaa)) - 4.45, nonzeros(dE4) .* data(:, 2), '.-', 'DisplayName', 'quartic fir for eff. V(S) - ED')
    %plot(abs(data(:,1)), data(:, 7) .* data(:, 2), 'o-', 'DisplayName', '1D Milnikov w/ Int')
    plot(-M_data(:, 1) - 4.45, M_data(:, 3) - M_data(:, 2), '*-', 'DisplayName', '3 particle ED')
    legend
    grid on
    %set(gca, 'Yscale', 'log')
    hold off

    %%
    figure(7)
    clf(figure(7))
    hold on
    beta = nonzeros(new_beta);
    alpha = nonzeros(new_alpha);
yline(0)
xline(-4.5)
    alphaa = nonzeros(alphaa)
    
    plot(nonzeros(alphaa), nonzeros(new_beta), 'DisplayName', 'const')
    plot(nonzeros(alphaa), nonzeros(new_alpha), 'DisplayName', 'alpha')

legend

    hold off