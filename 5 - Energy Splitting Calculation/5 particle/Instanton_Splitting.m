clc
clear all

disp('Tunneling Splitting calculation')

addpath('Data')
addpath('')
FS = 1;

for ind = 3:17
    
    name_ind = ind;
    name_str = ['TrajectoryData_5p_' num2str(name_ind)];

    load_data = load(name_str);
    data = load_data.IterData;
    % Number of Particles & Physicalvariables -----------------------------
    NoP     = data.NumberOfParticles;
    Alpha   = data.AlphaValues(ind);
    Eta     = data.CoulombStrength;
    Nt      = data.NumberOfPoinsInTrajectory;
    R       = data.RValues(ind);

    disp(['\alpha = ' num2str(Alpha)])


    % Euqilibrium positions -----------------------------------------------
    EquilibriumPositions = data.EquilibriumPositions(:, ind);

    %Trajectory -----------------------------------------------------------
    Trajectory = data.Trajectories;
    for i = 1:5
        Trajectory(i, :) = smooth(Trajectory(i, :));
    end
    
    % Time ----------------------------------------------------------------
    z   = data.time;
    dz  = data.dt;
    
    if FS
        figure(1)
        clf(figure(1))
        hold on
        for i = 1:NoP
            Trajectory(i, :) = smooth(Trajectory(i, :));
            plot(z, Trajectory(i, :), 'k.-')
            %plot(z, flip(Trajectory(i, :)), 'b.-')
        end
        yline(EquilibriumPositions, 'r')
        yline(-EquilibriumPositions, 'r')
        xline(0, 'b')
        yline(0, 'b')
        title(['Trajectories w/ \alpha = ' num2str(Alpha) ' and w/ r value = ' num2str(R)])
        xlabel('z - time')
        ylabel('\chi')
        hold off
    end
    
    % Calculating the arc length paramterization:
    [S, VS] = ArcLengthParametrization(Trajectory, NoP, Nt, Alpha, Eta);
    shiftAL = min(VS);
    if FS
        figure(2)
        clf(figure(2))
        hold on
        plot(z, S)
        xline(0)
        title('S(z)')
        xlabel('z - time')
        ylabel('S')
        hold off
    end 
%
    if FS
        figure(3)
        clf(figure(3))
        hold on
        plot(S, VS - shiftAL, '.-')
        xline(0)
        title('Arc length parametrized effective potential')
        xlabel('Arc length')
        ylabel('V(s)')
        hold off
    end 

    % Calculate the "velocity" of the trajectory
    [Velocity, VelocityUnitVector] = TrajectoryDifferentiation(Trajectory, NoP, Nt, z, dz, 1);

    if FS
        figure(4)
        clf(figure(4))
        hold on
        for i = 1:NoP
            plot(z, Velocity(i, :), 'k.-')
            plot(z, flip(Velocity(i, :)), 'b.-')
        end
        xline(0)
        title(['Velocities w/ \alpha = ' num2str(Alpha) ' and w/ r value = ' num2str(R)])
        xlabel('z - time')
        ylabel('\chi')
        hold off
    end
    %
    [fVector, DeltaVelocityUnitVector, DeltaPhi, RotationMatrix] = TrajectoryRotationMatrix(VelocityUnitVector, NoP, Nt);
    
    if FS
        figure(5)
        clf(figure(5))
        hold on
        for i = 1:Nt
            scalarproduct(i) = fVector(:, i)' * VelocityUnitVector(:, i);
        end
        plot(z, scalarproduct, '.-')
        xline(0)
        title('Scalar product of f and e')
        hold off
    end

    if FS
        figure(6)
        clf(figure(6))
        hold on
        plot(z, DeltaPhi)
        title('\Delta \phi (z)')
        hold off
    end
    %
    % Initializing tau vectors asa function of 
    [Tau0, EigVals, Tauspace] = TauspaceInit(Trajectory, Alpha, Velocity, Eta, NoP, Nt, RotationMatrix);
    %Tau0 = - Tau0;
    % Calculating the curvateure
    C = Curvature(Tauspace, z, dz, Velocity, NoP, Nt);

    if FS
        figure(7)
        clf(figure(7))
        hold on
        for i = 1:NoP
            plot(z,(C(i, :)))
        end
        title('Curvature')
        hold off
    end

    % B matrix:
    BMatrix = BSpringMatrixInit(Trajectory, Alpha, Tauspace, Eta, NoP);

    % Determine the Omega square matrix
    OmegaSquare = OmegaSquared(BMatrix, C, NoP, Nt);

    for i = 1:4
        for j = 1:4
            if i == j
                OmegaSquare(i, j, :) = (OmegaSquare(i, j, :));
            else
                %OmegaSquare(i, j, :) = 0;
            end
        end
    end

    if FS
        figure(8)
        clf(figure(8))
        hold on
        for i = 1:NoP-1
            for j = 1:NoP-1
                if i == j
                    plot(z, reshape(OmegaSquare(i, j, :), [Nt 1]), 'k.-')
                else
                    plot(z, reshape(OmegaSquare(i, j, :), [Nt 1]), 'b.-')
                end
            end
        end
        xline(0)
        title('Omega matrix')
        hold off
    end

    % Using the previous omega square matrix to solve the differential
    % equation:
    Xi = DifferentialEquationSolver(z, dz, NoP, Nt, OmegaSquare, R);
    
    if FS
        figure(9)
        clf(figure(9))
        hold on
        for i = 1:NoP-1
            for j = 1:NoP-1
                plot(z, reshape(Xi(i, j, :), [Nt 1]))
            end
        end
        xline(0)
        title('Xi_matrix elements')
        ylim([0 10])
        hold off
    end

    %Calculating everything:
    [PerpPropagator, Trace1, Trace2, Integ1, Integ2, Integ3] = Prefactor(Xi, OmegaSquare, EigVals, z, dz, Alpha, R, NoP, Nt);

    % Putting together everything:
    [shift, Spectra, SplitAL, Psi, VS, S] = Splitting(Xi, EigVals, z, Alpha, R, NoP, Nt, Eta, Trajectory, data.EnergySeries(end), S, VS);

    if FS
        figure(10)
        clf(figure(10))
        hold on
        title(SplitAL)
        plot(Spectra)
        %title('Curvature')
        hold off
    end
    SPLITTINGS(ind, 1) = SplitAL * PerpPropagator;
    SPLITTINGS(ind, 2) = Alpha;
    SPLITTINGS(ind, 3) = SplitAL;
end

figure(11)
clf(figure(11))
hold on
plot(SPLITTINGS(:, 2), SPLITTINGS(:, 1))
set(gca, 'Yscale', 'log')
hold off