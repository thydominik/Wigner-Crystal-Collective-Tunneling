clc
clear all

data = zeros(15,2);

for ind = 1:5:71
    Sf = [];
    % Pulling in the trajectory
    name = ind;
    nameSTR = ['P_' num2str(name)];

    %loading a trajectory that is previoulsy determined by a MC simulatiton
    trajectory_load              = load(nameSTR);%load('200_point_ztime_r_1_3_a_8_eta_20');

    %this must be the same eq_pos file that used in the MC, but therefore it's
    %not important bc the trajectory's endpoints will be the eq_positions!
    equilibrium_positions       = load('EqPos_eta20_alpha_5_20'); %4th column is the alpha value!
    trajectory                  = trajectory_load.position;

    %this will select the desired endpoints and alpha value
    state                       = ind;  
    eq_pos                      = equilibrium_positions.eqpos(:, state);

    %it's useful not to set in stone the division and particle number for later
    %more-particle cases
    [particle_n, N_division]    = size(trajectory);


    %This needs to be carried over from the trajectory calculation
    R = linspace(1.3, 0.3, 151);
    

    eps             = 10^-15;   %this have to be changed manually if it changes in the trajectory code!!!
    r               = R(state); %match this with the M.C. simulation
    alpha           = eq_pos(4); disp(['Alpha= ', num2str(alpha)])  %from eq_pos file!!!!
    eta             = 20;       %don't try to change this. Fitted value(experimetnts) = 18.813
    limits          = 50;       %1.35;  % +/- T

    z_time          = linspace(-1, 1, N_division);

    [chiS, S, VS]   = f_arclength(trajectory, alpha, eta, z_time);
    
    interpdiv = 0.005;
    dS = interpdiv;
    Sq = min(S):interpdiv:max(S);
    LS = length(Sq);
    S = S - max(S)/2;
    Sq= Sq - max(Sq)/2;
    VS_int = interp1(S, VS, Sq, 'Spline');
    figure(2)
    clf(figure(2))
    hold on
    plot(S, VS - min(VS), '.')
    plot(Sq, VS_int - min(VS_int), '.')
    hold off
    
    

    [gof, fc]       = f_fitting_VS(S(1:10), VS(1:10));
    extra = 0.5;
    for i = 1:(extra * LS)
        Sf(i) = (-i) * dS - max(Sq);
    end
    Sf = [flip(Sf) Sq];
    for i = 1:(extra * LS)
        Sf(end + 1) = i * dS + max(Sq);
    end

    figure(3)
    clf(figure(3))
    hold on
    plot(Sq)
    plot(Sf)
    hold off

    VSf = [(fc.b * (Sf(1:LS*extra).^2 - min(S)^2).^2) (VS_int - min(VS_int)) (fc.b * (Sf(end - LS*extra + 1:end).^2 - min(S)^2).^2)];
    
    figure(3)
    clf(figure(3))
    hold on
    plot(Sf, VSf)
    plot(S, VS - min(VS))
    hold off
    omegaS_sq_S     = 4 * fc.b * (3 * max(S)^2 - max(S)^2);

    [Spectra1] = Schrodinger_VS(VS, S, fc.a, fc.b, min(S)^2);
    [Psi2, Spectra2] = Schrodinger_VSF(Sf,fc.a + VSf);

    figure(4)
    clf(figure(4))
    hold on
    plot(Spectra1)
    plot(Spectra2)
    hold off
    dE1(ind) = Spectra1(2) - Spectra1(1);
    dE2(ind) = Spectra2(2) - Spectra2(1);
    a(ind) = alpha;
end
%%
M_data1 = load('E_Schrodinger_3e_eta_20.00_beta_0.01_N_100.dat');
data = load('WorkSpace.mat');
data = data.data;
figure(5)
clf(figure(5))
hold on
plot(-nonzeros(a), nonzeros(dE1) .* data(:, 2), '.-')
plot(-nonzeros(a), nonzeros(dE2) .* data(:, 2), 'o-')
plot(-M_data1(:, 1), M_data1(:, 3) - M_data1(:, 2), '.-', 'DisplayName', 'Beta = 0.01')
set(gca, 'Yscale', 'log')
hold off