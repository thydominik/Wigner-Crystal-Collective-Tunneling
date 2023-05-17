%% Supplementary material figures
clc
clear all
close all

set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultFigureWindowStyle','normal')

%% Figure 1: Tunneling trajectories

clear all
close all
clc

FontSize = 20;
Position = [1 1 14 10];

Figure = figure('color', 'white', 'units', 'inches', 'Position', Position);

hold on
% Fontsize:
    set(gca, 'FontSize', FontSize)

% Axis limits:
    %xlim([-5.5 5.5])
    %xticks([-5 -2.5 0 2.5 5])
    %ylim([-0.1 0.1])
    %yticks([-0.10 -0.05 0.0 0.05 0.10])
    %zlim([])
    %axis square

% Axis Labels:
    xlabel('$z$', 'Interpreter', 'latex', 'FontSize', FontSize)
    ylabel('$\chi(z)$', 'FontWeight', 'Bold', 'Interpreter', 'latex', 'FontWeight', 'Bold', 'FontSize', FontSize)
    %zlabel(sprintf('###'), 'Interpreter', 'latex', 'FontSize', FontSize)
% Data: 
    %addpath('D:\BME PhD\.Wigner Crystal Collective Tunneling\CollectiveTunneling\6 - Polarization\3 - Particle\WaveFunction data')
    Polarization = load('Polarization3Particle.mat'); Polarization = Polarization.Polarization
    for alphaInd = 1:26
        for kappaInd = 1:201
            name1 = ['X_' num2str(alphaInd) '_' num2str(kappaInd)];
            name2 = ['Y_' num2str(alphaInd) '_' num2str(kappaInd)];
            name3 = ['Z_' num2str(alphaInd) '_' num2str(kappaInd)];
            PsiX = load(name1); PsiX = PsiX.PsiX;
            PsiY = load(name2); PsiY = PsiY.PsiY;
            PsiZ = load(name3); PsiZ = PsiZ.PsiZ;
            LD1 = PsiX .* PsiX;
            LD2 = PsiY .* PsiY;
            LD3 = PsiZ .* PsiZ;
            LDOS(alphaInd, kappaInd, :) = LD1 + LD2 + LD3;
        end
    end
% Plotting:
    alphaIndicies = [11 14 17 18 19 20];
    k = 1;
    clim([0 0.9]);
    for i = alphaIndicies
        subplot(2, 3, k)
            alphaInd = i;
            hold on
            surf(linspace(-7, 7, 100), Polarization.Kappa, reshape(LDOS(alphaInd, :, :), [201 100]), 'EdgeColor', 'None', 'FaceColor', 'interp')
            plot3(Polarization.Polarization(alphaInd, :), Polarization.Kappa, 10 * ones(26,201),'k--', 'LineWidth', 3)
            xlim([-5 5])
            ax = gca;
            ax.FontSize = 14;
            %axis square
            title(['\alpha = ' num2str(Polarization.Alpha(alphaInd))], 'FontSize', 20)
            xlabel('\chi', 'FontSize', 20 + 5)
            ylabel('\epsilon', 'FontSize', 20 + 5)
            view([00 90])
            colormap turbo
            %colorbar
            lims = clim;
            clim([0 0.17]);
            hold off
            k = k + 1;
    end
% Colors:
    colormap turbo;
    c                   = colorbar('eastoutside'); %, 'Direction','reverse');
    c.Position = c.Position + [.05 .22 .01 .01];
    %c.Label.String      = '$\left| \Psi  \right|^2$';
    %c.Label.Interpreter = 'latex';
    %c.Label.FontSize    = FontSize + 5; 
% Additional Labels:
    %text( -5, 0.085, 10, '(c)', 'Color', 'white', 'interpreter', 'Latex', 'FontSize', FontSize);
% PaperSize:
    set(gcf,'paperunits','in');
    set(gcf,'papersize',[Position(3) + 1, Position(4) + 1]);
% Saving:
    fname   = sprintf('SupMatFig_DensitySeries.pdf');
    hfig    = gcf;
    print(hfig,'-bestfit','-dpdf', '-r960', fname);
%% Figure 2: Tunneling trajectories

clear all
close all
clc

FontSize = 18;
Position = [1 1 7 5];

Figure = figure('color', 'white', 'units', 'inches', 'Position', Position);

hold on
% Fontsize:
    set(gca, 'FontSize', FontSize)

% Axis limits:
    %xlim([-5.5 5.5])
    %xticks([-5 -2.5 0 2.5 5])
    %ylim([-0.1 0.1])
    yticks([-5 -2.5 0.0 2.5 5])
    %zlim([])
    %axis square
    box on
% Axis Labels:
    xlabel('$\vartheta$', 'Interpreter', 'latex', 'FontSize', FontSize + 5)
    ylabel('$\chi_i(\vartheta)$', 'FontWeight', 'Bold', 'Interpreter', 'latex', 'FontWeight', 'Bold', 'FontSize', FontSize + 5)
    %zlabel(sprintf('###'), 'Interpreter', 'latex', 'FontSize', FontSize)
% Data: 
    %addpath('D:\BME PhD\.Wigner Crystal Collective Tunneling\CollectiveTunneling\4 - Trajectory Calculation\3 particle\Standard MC Data')

    ThreeParticleTraj = load('Traj_3p_STDMC15.mat'); ThreeParticleTraj = ThreeParticleTraj.IterData
% Plotting:
    plot(linspace(-1, 1, 200), ThreeParticleTraj.Trajectories(2, :),'b-', 'LineWidth', 4.5)
    plot(linspace(-1, 1, 200), ThreeParticleTraj.Trajectories(1, :),'b-', 'LineWidth', 4.5)
    plot(linspace(-1, 1, 200), ThreeParticleTraj.Trajectories(3, :),'b-', 'LineWidth', 4.5)
    yline([ThreeParticleTraj.Trajectories(1, 1)-3*10^-2 ThreeParticleTraj.Trajectories(1, 200)+3*10^-2], '--', 'LineWidth', 2)
    yline([ThreeParticleTraj.Trajectories(2, 1)-3*10^-2 ThreeParticleTraj.Trajectories(2, 200)+3*10^-2], '--', 'LineWidth', 2)
    yline([ThreeParticleTraj.Trajectories(3, 1)-3*10^-2 ThreeParticleTraj.Trajectories(3, 200)+3*10^-2], '--', 'LineWidth', 2)
    yline(0, 'Alpha', 1, 'LineWidth', 2)
    xline(0, 'Alpha', 1, 'LineWidth', 2)
% Colors:
    %colormap turbo;
    %c                   = colorbar('eastoutside'); %, 'Direction','reverse');
    %c.Label.String      = '$\left| \Psi  \right|^2$';
    %c.Label.Interpreter = 'latex';
    %c.Label.FontSize    = FontSize + 5; 
% Additional Labels:
    %text( -1.6, 4.6, 0, '(a)', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize);
    text( 0.25, 1.10, 0, '$\alpha = 12.5$', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize + 3);
% PaperSize:
    set(gcf,'paperunits','in');
    set(gcf,'papersize',[Position(3) + 0.2, Position(4) + 0.2]);
% Saving:
    fname   = sprintf('Fig_3Particle_Trajectory.pdf');
    hfig    = gcf;
    print(hfig,'-bestfit','-dpdf', '-r960', fname);
%% Figure 3: Tunneling trajectories FIVE PARTICLES

clear all
close all
clc

FontSize = 20;
Position = [1 1 7 5];

Figure = figure('color', 'white', 'units', 'inches', 'Position', Position);

hold on
% Fontsize:
    set(gca, 'FontSize', FontSize)

% Axis limits:
    %xlim([-5.5 5.5])
    %xticks([-5 -2.5 0 2.5 5])
    ylim([-5.5 5.5])
    yticks([-5 -2.5 0.0 2.5 5])
    %zlim([])
    %axis square

% Axis Labels:
    xlabel('$z$', 'Interpreter', 'latex', 'FontSize', FontSize)
    ylabel('$\chi_i(z)$', 'FontWeight', 'Bold', 'Interpreter', 'latex', 'FontWeight', 'Bold', 'FontSize', FontSize)
    %zlabel(sprintf('###'), 'Interpreter', 'latex', 'FontSize', FontSize)
% Data: 
    %addpath('D:\BME PhD\.Wigner Crystal Collective Tunneling\CollectiveTunneling\4 - Trajectory Calculation\5 particle')
    FiveParticleTraj = load('TrajectoryData_5p_4.mat'); FiveParticleTraj = FiveParticleTraj.IterData
% Plotting:
    plot(linspace(-1, 1, 160), FiveParticleTraj.Trajectories(2, :),'b-', 'LineWidth', 4.5)
    plot(linspace(-1, 1, 160), FiveParticleTraj.Trajectories(1, :),'b-', 'LineWidth', 4.5)
    plot(linspace(-1, 1, 160), FiveParticleTraj.Trajectories(3, :),'b-', 'LineWidth', 4.5)
    plot(linspace(-1, 1, 160), FiveParticleTraj.Trajectories(4, :),'b-', 'LineWidth', 4.5)
    plot(linspace(-1, 1, 160), FiveParticleTraj.Trajectories(5, :),'b-', 'LineWidth', 4.5)
    yline([FiveParticleTraj.Trajectories(1, 1)-3*10^-2 FiveParticleTraj.Trajectories(1, 160)+3*10^-2], '--', 'LineWidth', 2)
    yline([FiveParticleTraj.Trajectories(2, 1)-3*10^-2 FiveParticleTraj.Trajectories(2, 160)+3*10^-2], '--', 'LineWidth', 2)
    yline([FiveParticleTraj.Trajectories(3, 1)-3*10^-2 FiveParticleTraj.Trajectories(3, 160)+3*10^-2], '--', 'LineWidth', 2)
    yline([FiveParticleTraj.Trajectories(4, 1)-3*10^-2 FiveParticleTraj.Trajectories(4, 160)+3*10^-2], '--', 'LineWidth', 2)
    yline([FiveParticleTraj.Trajectories(5, 1)-3*10^-2 FiveParticleTraj.Trajectories(5, 160)+3*10^-2], '--', 'LineWidth', 2)
    yline(0, 'Alpha', 1, 'LineWidth', 2)
    xline(0, 'Alpha', 1, 'LineWidth', 2)
% Colors:
    %colormap turbo;
    %c                   = colorbar('eastoutside'); %, 'Direction','reverse');
    %c.Label.String      = '$\left| \Psi  \right|^2$';
    %c.Label.Interpreter = 'latex';
    %c.Label.FontSize    = FontSize + 5; 
% Additional Labels:
    text( -0.95, 5, 0, '(#)', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize);
    text( 0.5, 0.85, 0, '$\alpha = 10.5$', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize);
% PaperSize:
    set(gcf,'paperunits','in');
    set(gcf,'papersize',[Position(3) + 1, Position(4) + 1]);
% Saving:
    fname   = sprintf('SupMatFig_5Particle_Trajectory.pdf');
    hfig    = gcf;
    print(hfig,'-bestfit','-dpdf', '-r960', fname);

%% Figure 4: Frequencies

clear all
close all
clc


FontSize = 19;
Position = [1 1 7 5];

Figure = figure('color', 'white', 'units', 'inches', 'Position', Position);

hold on
% Fontsize:
    set(gca, 'FontSize', FontSize)

% Axis limits:
    %xlim([-5.5 5.5])
    %xticks([-5 -2.5 0 2.5 5])
    %ylim([-5.5 5.5])
    yticks([0.0 2.0 4.0 6.0 8.0])
    ytickformat('%.1f')
    %zlim([])
    %axis square

% Axis Labels:
    xlabel('$\alpha$', 'Interpreter', 'latex', 'FontSize', FontSize+3)
    ylabel('$\omega_i \left( \alpha \right)$', 'FontWeight', 'Bold', 'Interpreter', 'latex', 'FontWeight', 'Bold', 'FontSize', FontSize+3)
    %zlabel(sprintf('###'), 'Interpreter', 'latex', 'FontSize', FontSize)
% Data: 
    F = load('Freqs.mat'); F = F.Frequencies;
    A = load('AlphaVals.mat'); A = A.AlphaValues;
% Plotting:
    plot(A, sqrt(F), '-', 'LineWidth', 3.5)
    xline(4.443, 'LineWidth', 2)
    %grid
    box
% Colors:
    %colormap turbo;
    %c                   = colorbar('eastoutside'); %, 'Direction','reverse');
    %c.Label.String      = '$\left| \Psi  \right|^2$';
    %c.Label.Interpreter = 'latex';
    %c.Label.FontSize    = FontSize + 5; 
% Additional Labels:
    %text( -1.35, 5, 0, '(#)', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize);
    %text( 0.5, 0.85, 0, '$\alpha = 10.5$', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize);
% PaperSize:
    set(gcf,'paperunits','in');
    set(gcf,'papersize',[Position(3) + 1, Position(4) + 1]);
% Saving:
    fname   = sprintf('Fig_Freqs.pdf');
    hfig    = gcf;
    print(hfig,'-bestfit','-dpdf', '-r960', fname);

%% Figure 5: effetivePotential

clear all
close all
clc

FontSize = 30;
Position = [1 1 7 5];

Figure = figure('color', 'white', 'units', 'inches', 'Position', Position);

hold on
% Fontsize:
    set(gca, 'FontSize', FontSize)

% Axis limits:
    xlim([-4 4])
    %xticks([-5 -2.5 0 2.5 5])
    ylim([0 25])
    %yticks([-5 -2.5 0.0 2.5 5])
    %zlim([])
    %axis square

% Axis Labels:
    xlabel('$s$', 'Interpreter', 'latex', 'FontSize', FontSize)
    ylabel('$v(s)$', 'FontWeight', 'Bold', 'Interpreter', 'latex', 'FontWeight', 'Bold', 'FontSize', FontSize)
    %zlabel(sprintf('###'), 'Interpreter', 'latex', 'FontSize', FontSize)
% Data: 
    %addpath('D:\BME PhD\.Wigner Crystal Collective Tunneling\CollectiveTunneling\4 - Trajectory Calculation\3 particle\Standard MC Data')

    ThreeParticleTraj = load('Traj_3p_STDMC15.mat'); ThreeParticleTraj = ThreeParticleTraj.IterData

    [S, VS]         = ArcLengthParametrization(ThreeParticleTraj.Trajectories, 3, 200, 12.5, 20);    % Arc length paramterization
    VS              = VS - min(VS);     % Shifting the effective potential to 0
    NoPS            = 300;              % Points in the interpolation
    Sq              = linspace(min(S), max(S), NoPS);   % New arc length parameter
    dS              = Sq(2) - Sq(1);    % difference between arc length points
    LSq             = length(Sq);       % NoPoints in the new S
    VS_interpolate  = interp1(S, VS, Sq, 'Spline');     % Interpolating
    VS_interpolate  = VS_interpolate - min(VS_interpolate);     % Shifting the interpolation to 0 (should be already at 0 tho)

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
% Plotting:
box
    plot(SS, VS_interpolate,'k-', 'LineWidth', 5)
    plot(-2.83, 0, 'ro', 'MarkerSize', 13, 'MarkerFaceColor', 'r')
    plot(2.83, 0, 'ro', 'MarkerSize', 13, 'MarkerFaceColor', 'r')
% Colors:
    %colormap turbo;
    %c                   = colorbar('eastoutside'); %, 'Direction','reverse');
    %c.Label.String      = '$\left| \Psi  \right|^2$';
    %c.Label.Interpreter = 'latex';
    %c.Label.FontSize    = FontSize + 5; 
% Additional Labels:
    %text( -3.8, 25, 0, '(#)', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize);
    %text( 2, 20, 0, '$\alpha = 10.5$', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize);
% PaperSize:
    set(gcf,'paperunits','in');
    set(gcf,'papersize',[Position(3) + 1, Position(4) + 1]);
% Saving:
    fname   = sprintf('SupMatFig_EffectivePotential.pdf');
    hfig    = gcf;
    print(hfig,'-bestfit','-dpdf', '-r960', fname);