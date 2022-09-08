clc
clear all
format long

FS = 1;         % Figure Switch = 0 - Off; 1 - On

disp('Instanton Prefactor and Tunneling splitting calculation')

addpath('D:\BME PhD\.Wigner Crystal Collective Tunneling\CollectiveTunneling\5 - Energy Splitting Calculation\5 particle\Trajectories')

%This loops goes through the previously calculated trajectories:
for ind = 1:13
    %Clearing a variable in each cycle: This will hold some of the arc
    %length
    Sf = [];

    % Pulling in the trajectory
    name = ind;

    nameSTR = ['TrajectoryData_5P_' num2str(name)];

    %loading a trajectory that is previoulsy determined by a MC simulatiton
    trajectory_load = load(nameSTR);
    trajecotry_data = trajectory_load.IterData;

    % Separating the equilibrium positions
    equilibrium_positions   = trajecotry_data.EquilibriumPositions(:, ind);
    Trajectory              = trajecotry_data.Trajectories;

    PN  = 5;
    NoP = length(Trajectory(1, :));

    disp(['Particle # = ',num2str(PN)])
    disp(['Trajectory points = ', num2str(NoP)])

    % z imaginery time paramter
    if ind <= 5
        r       = trajecotry_data.RValues(ind);
        Alpha   = trajecotry_data.AlphaValues(ind);
    else
        r       = trajecotry_data.RValues(ind);
        Alpha   = trajecotry_data.AlphaValues(ind);
    end
    
    Eta     = 20;       % Change this if u r a Heretic!

    z = trajecotry_data.time;

    if FS
        figure(1)
        clf(figure(1))
        hold on
        title('Trajectories')
        xlabel('z time')
        ylabel('\chi_i')
        xline(0)
        yline(0)
        yline(equilibrium_positions)
        yline(-equilibrium_positions)
        plot(z, Trajectory(1, :), 'o-')
        plot(z, Trajectory(2, :), 'o-')
        plot(z, Trajectory(3, :), 'o-')
        plot(z, Trajectory(4, :), 'o-')
        plot(z, Trajectory(5, :), 'o-')
        hold off
    end

    % Using the trajectories, creating the arc length param and the arc
    % length parametrized V(S) potential.
    [S, VS]         = ArcLengthParametrization(Trajectory, PN, NoP, Alpha, Eta);    % Arc length paramterization
    VS              = VS - min(VS);     % Shifting the effective potential to 0
    NoPS            = 300;              % Points in the interpolation
    Sq              = linspace(min(S), max(S), NoPS);   % New arc length parameter
    dS              = Sq(2) - Sq(1);    % difference between arc length points
    LSq             = length(Sq);       % NoPoints in the new S
    VS_interpolate  = interp1(S, VS, Sq, 'Spline');     % Interpolating
    VS_interpolate  = VS_interpolate - min(VS_interpolate);     % Shifting the interpolation to 0 (should be already at 0 tho)

    if FS
        figure(2)
        clf(figure(2))
        hold on
        title('Effective potential and its interpolation')
        plot(S, VS - min(VS), 'o-')
        plot(Sq, VS_interpolate, '-', 'LineWidth', 3)
        hold off
    end

    % Calculating the Arc length 1D Schrödinger problem
    % Fitting the V(S) potential with a quartic potential
    [gof, fc] = f_fitting_VS_2(S(1:20), VS(1:20)); % Fit the first part of the potential

    NoNP = 1000; %The Number of New Points in V(S) (to continue the potential)
    SS = Sq;
    for i = 1:NoNP
        VS_interpolate = [fc.b * ((-i*dS + min(Sq)) - min(Sq))^2 VS_interpolate];
        SS = [(-i*dS + min(Sq)) SS];
        VS_interpolate = [VS_interpolate fc.b * ((i*dS + max(Sq)) - max(Sq))^2];
        SS = [SS (i*dS + max(Sq))];
    end
    SS = SS - min(SS);
    SS = SS - max(SS)/2;

    % --------------------------------
    NoPS = NoPS + 2*NoNP;
    Hami        = zeros(NoPS, NoPS);
    Kinetic     = zeros(NoPS, NoPS);
    Potential   = zeros(NoPS, NoPS);

    K = 1/(2 * dS^2);

    for i = 2:NoPS
        Kinetic(i - 1, i) = -K;
        Kinetic(i, i - 1) = -K;
    end

    for i = 1:NoPS
        Potential(i, i) = VS_interpolate(i);
    end

    if FS
        figure(3)
        clf(figure(3))
        hold on
        title('Interpolation of the Effective potential (and the continuation of it)')
        plot(SS, VS_interpolate, 'k.')
        plot(S, VS, 'O')
        hold off
    end

    Hamilton        = Potential + Kinetic;
    [Psi, Spectra]  = eig(Hamilton);
    Spectra         = diag(Spectra);
    SplitED(ind)    = Spectra(2) - Spectra(1);
    if FS
        figure(4)
        clf(figure(4))
        hold on
        title(["Wavefunctions and the splitting: " num2str(Spectra(2) - Spectra(1))])
        plot(SS, VS_interpolate/norm(VS_interpolate))
        plot(SS, Psi(:, 1), '.-')
        plot(SS, Psi(:, 2), 'r.-')
        yline(0)
        xline(0)
        hold off
    end


    %In order to get the tangent vector first calc the derivative of the curve
    dz      = z(2) - z(1);
    [Velocity, VelocityUnitVector] = TrajectoryDifferentiation(Trajectory, PN, NoP, z, dz, 0);

    if FS
        figure(5)
        clf(figure(5))
        hold on
        title('Derivative of the Trajectory')
        title(['Velocities w/ \alpha = ' num2str(Alpha) ' and w/ r value = ' num2str(r)])
        plot(z, Velocity)
        xlabel('z')
        ylabel('d/dz \chi')
        hold off
    end

    if FS
        figure(6)
        clf(figure(6))
        hold on
        title('Velocity unit vector components')
        plot(z, VelocityUnitVector)
        xline(0)
        hold off
    end

    [fVector, DeltaVelocityUnitVector, DeltaPhi, RotationMatrix] = TrajectoryRotationMatrix(VelocityUnitVector, PN, NoP);
    if FS
        figure(7)
        clf(figure(7))
        hold on
        scalarproduct = [];
        for i = 1:NoP-1
            scalarproduct(i) = fVector(:, i)' * VelocityUnitVector(:, i);
        end
        plot(z(1:end-1), scalarproduct, '.-')
        xline(0)
        title('Scalar product of f and e')
        hold off
    end

    if FS
        figure(8)
        clf(figure(8))
        hold on
        plot(z, DeltaPhi)
        title('\Delta \phi (z)')
        hold off
    end

    % Initializing tau vectors asa function of
    [Tau0, EigVals, Tauspace] = TauspaceInit(Trajectory, Alpha, Velocity, Eta, PN, NoP, RotationMatrix);
    %Tau0 = - Tau0;
    % Calculating the curvateure
    C = Curvature(Tauspace, z, dz, Velocity, PN, NoP);

    if FS
        figure(9)
        clf(figure(9))
        hold on
        for i = 1:PN
            plot(z,(C(i, :)))
        end
        title('Curvature')
        hold off
    end

    % B matrix:
    BMatrix = BSpringMatrixInit(Trajectory, Alpha, Tauspace, Eta, PN);

    % Determine the Omega square matrix
    OmegaSquare = OmegaSquared(BMatrix, C, PN, NoP);

    for i = 1:PN-1
        for j = 1:PN-1
            if i == j
                OmegaSquare(i, j, :) = (OmegaSquare(i, j, :));
            else
                %OmegaSquare(i, j, :) = 0;
            end
        end
    end

    if FS
        figure(10)
        clf(figure(10))
        hold on
        for i = 1:PN-1
            for j = 1:PN-1
                if i == j
                    plot(z, reshape(OmegaSquare(i, j, :), [NoP 1]), 'k.-')
                else
                    plot(z, reshape(OmegaSquare(i, j, :), [NoP 1]), 'b.-')
                end
            end
        end
        xline(0)
        title('Omega matrix')
        hold off
    end

    % Using the previous omega square matrix to solve the differential
    % equation:
    Xi = DifferentialEquationSolver(z, dz, PN, NoP, OmegaSquare, r);

    if FS
        figure(11)
        clf(figure(11))
        hold on
        for i = 1:PN-1
            for j = 1:PN-1
                plot(z, reshape(Xi(i, j, :), [NoP 1]))
            end
        end
        xline(0)
        title('Xi_matrix elements')
        ylim([0 10])
        hold off
    end

    %Calculating everything:
    [PerpPropagator, Trace1, Trace2, Integ1, Integ2, Integ3] = Prefactor(Xi, OmegaSquare, EigVals, z, dz, Alpha, r, PN, NoP);

    % Putting together everything:
    [Split_EasyMilnikov, Split_Landau, Split_Instanton, Split_Prop_a, Split_Prop_b] = Splitting(Xi, EigVals, z, Alpha, r, PN, NoP, Eta, Trajectory, trajecotry_data.Action, S, VS);
    
    if FS
        figure(12)
        clf(figure(12))
        hold on
        plot(Spectra)
        %title('Curvature')
        hold off
    end

    SPLITTINGS(ind, 1) = Alpha;
    SPLITTINGS(ind, 2) = SplitED(end);
    SPLITTINGS(ind, 3) = SplitED(end) * PerpPropagator;
    SPLITTINGS(ind, 4) = Split_EasyMilnikov * PerpPropagator;
    SPLITTINGS(ind, 5) = Split_Landau * PerpPropagator;
    SPLITTINGS(ind, 6) = Split_Instanton * PerpPropagator;
    SPLITTINGS(ind, 7) = Split_Prop_a * PerpPropagator;
    SPLITTINGS(ind, 8) = Split_Prop_b * PerpPropagator;
    SPLITTINGS(ind, 9) = PerpPropagator;
    
end

save('Standard5particleSplitting', 'SPLITTINGS')

%%
% SIMPMC = load('SimplifiedSplitting.mat');
STDMC = load('Standard5particleSplitting.mat');
% RESTMC = load('RestrictedSplitting.mat');
% 
ED = load('EDSplitting_5_particles_restricted_Nx1_15_Nx2_20_Nx3_40_Nx4_20_Nx5_15_beta_1e-05.mat');
ED = ED.data;
ED_a = ED(:, 1);
ED_d = ED(:, 2);
% 
% DMRG = load('Delta_E_DMRG_Norb_8_eta_20.00.mat');
% DMRG_a = DMRG.alpha_list;
% DMRG_d = DMRG.delta_E_DMRG;

figure(13)
clf(figure(13))
hold on
xlabel('\alpha', 'FontSize', 20)
ylabel('\Delta', 'FontSize', 20)
title('5 Particle Splitting')
legend
%plot(SPLITTINGS(:, 1), SPLITTINGS(:, 2), 'o-', 'DisplayName', '1D part from ED')
plot(STDMC.SPLITTINGS(:, 1), STDMC.SPLITTINGS(:, 3), '.-', 'LineWidth', 2, 'DisplayName', '1D part ED * N-1 Milnikov - STD')
%plot(RESTMC.SPLITTINGS(:, 1), RESTMC.SPLITTINGS(:, 3), '.-', 'LineWidth', 2, 'DisplayName', '1D part ED * N-1 Milnikov - REST')
%plot(SIMPMC.SPLITTINGS(:, 1), SIMPMC.SPLITTINGS(:, 3), '.-', 'LineWidth', 2, 'DisplayName', '1D part ED * N-1 Milnikov - SIMP')
%plot(SPLITTINGS(:, 1), SPLITTINGS(:, 4), '.-', 'DisplayName', 'Easy Milnikov')
%plot(SPLITTINGS(:, 1), SPLITTINGS(:, 5), '.-', 'DisplayName', 'Landau')
%plot(SPLITTINGS(:, 1), SPLITTINGS(:, 6), '.-', 'DisplayName', 'Coleman')
%plot(SPLITTINGS(:, 1), SPLITTINGS(:, 7), 'o-', 'DisplayName', 'Integral ver.1')
%plot(SPLITTINGS(:, 1), SPLITTINGS(:, 8), 'o-', 'DisplayName', 'Integral ver. 2')
plot(ED(:, 1), ED(:, 2), 's-', 'DisplayName', 'ED')
%plot(-DMRG_a, DMRG_d, '--', 'DisplayName', 'DMRG')

%set(gca, 'Yscale', 'log')
xlim([7.8 14])
ylim([10^-5 3])
xline(10.857299871771385, 'r', 'DisplayName', 'Elbow point')
xline(7.84, 'b', 'DisplayName', 'Critical Alpha')
%xline(4.45, 'k', 'DisplayName', 'Classical critical value')
hold off
    disp('---------------------------------000-------------------------------------')



