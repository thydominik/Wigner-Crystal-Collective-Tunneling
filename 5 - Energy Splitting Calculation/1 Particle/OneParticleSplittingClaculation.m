clc
clear all
format long

%% Loading in the calculated trajecotires:

addpath('D:\BME PhD\.Wigner Crystal Collective Tunneling\CollectiveTunneling\4 - Trajectory Calculation\1 particle\Data')

for stateInd = 1:19
    NameSTR     = ['Traj_1p_' num2str(stateInd)];
    TrajLoad    = load(NameSTR);
    TrajState = TrajLoad.IterData;

    Action(stateInd)        = TrajState.Action;
    z                       = TrajState.time;
    Trajectory(stateInd, :) = TrajState.Trajectories;
    Alpha(stateInd)         = TrajState.AlphaValues(stateInd);
    R(stateInd)             = TrajState.RValue;
end

%% Calculating the frequency

for stateInd = 1:19
    Omega(stateInd) = sqrt(2 * abs(Alpha(stateInd)));
end

%%

for i = 1:19

    %Potential:
    V = 0.25 .* (Trajectory(i, :).^2 - abs(Alpha(i))).^2;

    %Classical Impulse:
    P   = sqrt(2 * V);
    P_0 = sqrt(2 * 0.25 * Alpha(i)^2);

    % (21) integral
    B_int = 0;

    B_pref = sqrt(4 * Omega(i) / pi) * P_0;

   %Derivative of P:
    dP(1) = (P(2) - P(1))./(Trajectory(i, 2) - Trajectory(i, 1));
    for k = 2:length(P)-1
        dP(k) = (P(k+1) - P(k-1))./(Trajectory(i, k + 1) - Trajectory(i, k - 1));
    end
    dP(length(P)) = (P(end) - P(end-1))./(Trajectory(i, end) - Trajectory(i, end-1));


    % Integrand
    temp_func = (R(i)./(1 - z.^2)) .* (Omega(i) - dP);
    for intInd = 2:100
        B_int = B_int + (z(intInd+1) - z(intInd)) * 0.5 * (temp_func(intInd + 1) + temp_func(intInd))
    end
    
    NumericalSplitting(i)   =  B_pref * exp(-Action(i)) * exp(B_int);
    Action_a                = 2/3 * sqrt(2) * Alpha(i).^(3/2);
    AnalyiticalSplitting(i) = sqrt(Omega(i)/pi) * 4 * Omega(i) * sqrt(Alpha(i)) * exp(-Action_a);
end

%%
DMRGdata = load("OneParticleDMRG.txt");
EDdata      = load('EDSplitting_1_particle_Nx_800.mat');

figure(1)
clf(figure(1))
hold on
plot(Alpha, NumericalSplitting, 'o-', 'DisplayName', 'Numerical Calculation (from trajectories)')
plot(Alpha, AnalyiticalSplitting, 'd-', 'DisplayName', 'Analytical calculation')
plot(-DMRGdata(:, 1), DMRGdata(:, 2), '*-', 'DisplayName', 'DMRG')
plot(EDdata.data(1:5:end, 1), EDdata.data(1:5:end, 2), 's-', 'DisplayName', 'ED')
%set(gca, "Yscale", "log")
xline(2.934, 'DisplayName',"Elbow's end")
grid
legend
xlabel('\alpha', 'FontSize', 20)
ylabel('\Delta', 'FontSize', 20)
hold off

figure(2)
clf(figure(2))
hold on
title('Numerical calculation and ED ratio')
xlabel('\alpha', 'FontSize', 20)
ylabel('\Delta diff.', 'FontSize', 20)
grid
legend
plot(Alpha, abs(EDdata.data(:, 1)./ AnalyiticalSplitting'), 'o-', 'DisplayName', 'ED/Num.')
set(gca, "Yscale", "log")
xline(2.934, 'DisplayName',"Elbow's end")
hold off