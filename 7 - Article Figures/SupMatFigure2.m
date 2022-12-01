clc
clear all

addpath('D:\BME PhD\.Wigner Crystal Collective Tunneling\CollectiveTunneling\4 - Trajectory Calculation\3 particle\Standard MC Data')

ThreeParticleTraj = load('Traj_3p_STDMC15.mat'); ThreeParticleTraj = ThreeParticleTraj.IterData

figure(12)
clf(figure(12))
hold on
plot(linspace(-1, 1, 200), ThreeParticleTraj.Trajectories(2, :),'b-', 'LineWidth', 2)
plot(linspace(-1, 1, 200), ThreeParticleTraj.Trajectories(1, :),'b-', 'LineWidth', 2)
plot(linspace(-1, 1, 200), ThreeParticleTraj.Trajectories(3, :),'b-', 'LineWidth', 2)
yline([ThreeParticleTraj.Trajectories(1, 1) ThreeParticleTraj.Trajectories(1, 200)])
yline([ThreeParticleTraj.Trajectories(2, 1) ThreeParticleTraj.Trajectories(2, 200)])
yline([ThreeParticleTraj.Trajectories(3, 1) ThreeParticleTraj.Trajectories(3, 200)])
yline(0)
xline(0)
box
axis square
xlabel('z', 'FontSize', 20)
ylabel('\chi(z)', 'FontSize', 20)
title(['\alpha = ' num2str(ThreeParticleTraj.AlphaValues(15))], 'FontSize', 20)
hold off


addpath('D:\BME PhD\.Wigner Crystal Collective Tunneling\CollectiveTunneling\4 - Trajectory Calculation\5 particle')
FiveParticleTraj = load('TrajectoryData_5p_4.mat'); FiveParticleTraj = FiveParticleTraj.IterData
figure(13)
clf(figure(13))
hold on
hold on
plot(linspace(-1, 1, 160), FiveParticleTraj.Trajectories(2, :),'b-', 'LineWidth', 2)
plot(linspace(-1, 1, 160), FiveParticleTraj.Trajectories(1, :),'b-', 'LineWidth', 2)
plot(linspace(-1, 1, 160), FiveParticleTraj.Trajectories(3, :),'b-', 'LineWidth', 2)
plot(linspace(-1, 1, 160), FiveParticleTraj.Trajectories(4, :),'b-', 'LineWidth', 2)
plot(linspace(-1, 1, 160), FiveParticleTraj.Trajectories(5, :),'b-', 'LineWidth', 2)
yline([FiveParticleTraj.Trajectories(1, 1) FiveParticleTraj.Trajectories(1, 160)])
yline([FiveParticleTraj.Trajectories(2, 1) FiveParticleTraj.Trajectories(2, 160)])
yline([FiveParticleTraj.Trajectories(3, 1) FiveParticleTraj.Trajectories(3, 160)])
yline([FiveParticleTraj.Trajectories(4, 1) FiveParticleTraj.Trajectories(4, 160)])
yline([FiveParticleTraj.Trajectories(5, 1) FiveParticleTraj.Trajectories(5, 160)])
yline(0)
xline(0)
box
axis square
xlabel('z', 'FontSize', 20)
ylabel('\chi(z)', 'FontSize', 20)
title(['\alpha = ' num2str(FiveParticleTraj.AlphaValues(4))], 'FontSize', 20)
hold off

    [S, VS]         = ArcLengthParametrization(ThreeParticleTraj.Trajectories, 3, 200, 12.5, 20);    % Arc length paramterization
    VS              = VS - min(VS);     % Shifting the effective potential to 0
    NoPS            = 300;              % Points in the interpolation
    Sq              = linspace(min(S), max(S), NoPS);   % New arc length parameter
    dS              = Sq(2) - Sq(1);    % difference between arc length points
    LSq             = length(Sq);       % NoPoints in the new S
    VS_interpolate  = interp1(S, VS, Sq, 'Spline');     % Interpolating
    VS_interpolate  = VS_interpolate - min(VS_interpolate);     % Shifting the interpolation to 0 (should be already at 0 tho)

    if 1
        figure(20)
        clf(figure(20))
        hold on
        title('Effective potential and its interpolation')
        plot(S, VS - min(VS), 'o-')
        plot(Sq, VS_interpolate, '-', 'LineWidth', 3)
        hold off
    end

    % Calculating the Arc length 1D Schrödinger problem
    % Fitting the V(S) potential with a quartic potential
    [gof, fc] = f_fitting_VS_2(S(1:20), VS(1:20)); % Fit the first part of the potential

    NoNP = 50; %The Number of New Points in V(S) (to continue the potential)
    SS = Sq;
    for i = 1:NoNP
        VS_interpolate = [fc.b * ((-i*dS + min(Sq)) - min(Sq))^2 VS_interpolate];
        SS = [(-i*dS + min(Sq)) SS];
        VS_interpolate = [VS_interpolate fc.b * ((i*dS + max(Sq)) - max(Sq))^2];
        SS = [SS (i*dS + max(Sq))];
    end
    SS = SS - min(SS);
    SS = SS - max(SS)/2;
%%
figure(14)
clf(figure(14))
hold on
plot(SS, VS_interpolate,'k-', 'LineWidth', 3)
plot(-2.83, 0, 'ro', 'MarkerSize', 7, 'MarkerFaceColor', 'r')
plot(2.83, 0, 'ro', 'MarkerSize', 7, 'MarkerFaceColor', 'r')
xlabel('S', 'FontSize', 20)
ylabel('V(S)', 'FontSize', 20)
axis square
hold off