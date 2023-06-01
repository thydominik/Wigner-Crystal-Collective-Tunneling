%% clc
clc
clear all
close all

set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultFigureWindowStyle','normal')
addpath('./Experimental data/')
addpath('./DMRG data/')
addpath('./WaveFunction data/')
%% Figure 1: Polarization 2 dimensional plot
%
clear all
close all
clc

FontSize = 20;
Position = [1 1 6 5];

Figure = figure('color', 'white', 'units', 'inches', 'Position', Position);

hold on
% Fontsize:
set(gca, 'FontSize', FontSize)

% Axis limits:
xlim([-5.5 5.5])
xticks([-5 -2.5 0 2.5 5])
ylim([-0.1 0.1])
yticks([-0.10 -0.05 0.0 0.05 0.10])
%yticks([-0.10  0.0  0.10])
%zlim([])
%axis square

% Axis Labels:
text( 3.75, 0.085, 10, '$\rho(\chi)$', 'Color', 'white', 'interpreter', 'Latex', 'FontSize', FontSize);
title('$\rho(\chi)$', 'Interpreter', 'latex', 'FontSize', FontSize)
xlabel('$\chi$', 'Interpreter', 'latex', 'FontSize', FontSize)
ylabel('$\epsilon$', 'FontWeight', 'Bold', 'Interpreter', 'latex', 'FontWeight', 'Bold', 'FontSize', FontSize + 5)
%zlabel(sprintf('###'), 'Interpreter', 'latex', 'FontSize', FontSize)
% Data:
Pol = load('Polarization3Particle.mat'); Pol = Pol.Polarization;
% Plotting:
AlphaInd = 20;

for alphaInd = AlphaInd:AlphaInd
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

surf(linspace(-7, 7, 100), Pol.Kappa, reshape(LDOS(AlphaInd, :, :), [201 100]), 'EdgeColor', 'none', 'FaceColor', 'interp')
plot3(Pol.Polarization(AlphaInd, :), Pol.Kappa, 10 * ones(26,201), 'k--', 'LineWidth', 5)
% Colors:
colormap turbo;
c                   = colorbar('eastoutside'); %, 'Direction','reverse');
%c.Label.String      = '$\rho(\chi)$';
%c.Label.Interpreter = 'latex';
%c.Label.FontSize    = FontSize + 5;
% Additional Labels:
%text( -5, 0.085, 10, '(a)', 'Color', 'white', 'interpreter', 'Latex', 'FontSize', FontSize);
% PaperSize:
set(gcf,'paperunits','in');
set(gcf,'papersize',[Position(3) + 1, Position(4) + 1]);
% Saving:
fname   = sprintf('Fig_Polarization_2D.pdf');
hfig    = gcf;
print(hfig,'-bestfit','-dpdf', '-r960', fname);

%% Figure 2: Wavefunction in potential POSITIVE epsilon:
%
clear all
close all
clc

FontSize = 30;
Position = [1 1 6 4.5];

Figure = figure('color', 'white', 'units', 'inches', 'Position', Position);

hold on
% Fontsize:
set(gca, 'FontSize', FontSize)

% Axis limits:
xlim([-5.5*160 5.5*160])
xticks([-500 0 500])
ylim([-0.001 0.003])
yticks([0 0.5])
%zlim([])
%axis square

% Axis Labels:
xlabel('$z [{\rm{nm}}]$', 'Interpreter', 'latex', 'FontSize', FontSize)
label = ylabel('$\rho(z) \left[ \frac{1}{nm}\right] $', 'FontWeight', 'Bold', 'Interpreter', 'latex', 'FontWeight', 'Bold', 'FontSize', FontSize + 5, 'Color','r')
%label.Position(1) = -5.8;
%label.Position(2) = 70;
%set(label, 'Rotation', 0)
%zlabel(sprintf('###'), 'Interpreter', 'latex', 'FontSize', FontSize)
% Data:
x = linspace(-7, 7, 1000);

Pol = load('Polarization3Particle.mat'); Pol = Pol.Polarization;

AlphaInd = 20;

WF1 = load(['X_' num2str(AlphaInd) '_1']); WF1 = WF1.PsiX; WF1 = interp1(linspace(-7, 7, 100), WF1, linspace(-7, 7, 1000), 'spline');
WF2 = load(['Y_' num2str(AlphaInd) '_1']); WF2 = WF2.PsiY; WF2 = interp1(linspace(-7, 7, 100), WF2, linspace(-7, 7, 1000), 'spline');
WF3 = load(['Z_' num2str(AlphaInd) '_1']); WF3 = WF3.PsiZ; WF3 = interp1(linspace(-7, 7, 100), WF3, linspace(-7, 7, 1000), 'spline');
WaveFunction = WF1 + WF2 + WF3;
% Plotting:
mult = 10000;
% Potential:
LinePlot = plot(x * 160, 1/mult * ((Pol.Kappa(1) * x) + 0.012 * (x.^2 - Pol.Alpha(AlphaInd)).^2) / 0.03, 'k-', 'LineWidth', 5);
LinePlot.Color(4) = 0.5;

yline(0, 'k-', 'LineWidth', 2.5, 'Alpha', 1);

% Wavefunctions:
trshld = 8*10^-4;
shift1 = -9.5;

%plot(x(abs(WF1) > trshld), mult * abs(WF1(abs(WF1) > trshld)), 'r-', 'LineWidth', 2.0)
%plot(x(abs(WF2) > trshld), mult * abs(WF2(abs(WF2) > trshld)), 'r-', 'LineWidth', 2.0)
%plot(x(abs(WF3) > trshld), mult * abs(WF3(abs(WF3) > trshld)), 'r-', 'LineWidth', 2.0)

plot(x*160, WaveFunction/160, 'r-', 'LineWidth', 6)
% Colors:
%colormap turbo;
%c                   = colorbar('eastoutside'); %, 'Direction','reverse');
%c.Label.String      = '$\left| \Psi  \right|^2$';
%c.Label.Interpreter = 'latex';
%c.Label.FontSize    = FontSize + 5;
% Additional Labels:
% Note that in the matlab figure, the label (a) is out of the fram,
% after saving the figure as pdf the label will be in a good place i
% think.
%text( -8.7, 65, 0, '(b)', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize);
text( 0, 0.45, 0, '$\epsilon < 0$','HorizontalAlignment', 'center', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize + 15);
% PaperSize:
set(gcf,'paperunits','in');
set(gcf,'papersize',[Position(3) + 1, Position(4) + 1]);
% Saving:
fname   = sprintf('Fig_WaveFunction_Potential_PosEps.pdf');
hfig    = gcf;
print(hfig,'-bestfit','-dpdf', '-r960', fname);
%% Figure 3: Wavefunction in potential NEGATIVE epsilon:
%
clear all
close all
clc

Version = 1; % 1 or 2
FontSize = 30;
Position = [1 1 6 4.5];

Figure = figure('color', 'white', 'units', 'inches', 'Position', Position);

hold on
% Fontsize:
set(gca, 'FontSize', FontSize)

% Axis limits:
xlim([-5.5*160 5.5*160])
xticks([-500 0 500])
ylim([-.1 0.5])
yticks([0 0.5 60])
%zlim([])
%axis square

% Axis Labels:
xlabel('$z [{\rm{nm}}]$', 'Interpreter', 'latex', 'FontSize', FontSize)
ylabel('$\rho(z) \left[ \frac{1}{nm}\right] $', 'FontWeight', 'Bold', 'Interpreter', 'latex', 'FontWeight', 'Bold', 'FontSize', FontSize + 5, 'Color','r')
%zlabel(sprintf('###'), 'Interpreter', 'latex', 'FontSize', FontSize)
% Data:
x = linspace(-7, 7, 1000);

Pol = load('Polarization3Particle.mat'); Pol = Pol.Polarization;

AlphaInd = 20;

WF1 = load(['X_' num2str(AlphaInd) '_1']); WF1 = WF1.PsiX; WF1 = interp1(linspace(-7, 7, 100), WF1, linspace(-7, 7, 1000), 'spline');
WF2 = load(['Y_' num2str(AlphaInd) '_1']); WF2 = WF2.PsiY; WF2 = interp1(linspace(-7, 7, 100), WF2, linspace(-7, 7, 1000), 'spline');
WF3 = load(['Z_' num2str(AlphaInd) '_1']); WF3 = WF3.PsiZ; WF3 = interp1(linspace(-7, 7, 100), WF3, linspace(-7, 7, 1000), 'spline');
WaveFunction = WF1 + WF2 + WF3;
% Plotting:
mult = 100;
% Potential:
LinePlot = plot(x*160, 1/mult * ((-Pol.Kappa(1) * x) + 0.012 * (x.^2 - Pol.Alpha(AlphaInd)).^2) / 0.03, 'k-', 'LineWidth', 5);
LinePlot.Color(4) = 0.5;

yline(0, 'k-', 'LineWidth', 2.5, 'Alpha', 1);

% Wavefunctions:
trshld = 8*10^-4;
shift1 = -9.5;

%plot(x(abs(WF1) > trshld), flip(mult * abs(WF1(abs(WF1) > trshld))), 'r-', 'LineWidth', 2.0)
%plot(x(abs(WF2) > trshld), flip(mult * abs(WF2(abs(WF2) > trshld))), 'r-', 'LineWidth', 2.0)
%plot(x(abs(WF3) > trshld), flip(mult * abs(WF3(abs(WF3) > trshld))), 'r-', 'LineWidth', 2.0)
plot(x*160, flip(WaveFunction), 'r-', 'LineWidth', 6)
if Version == 2
    y1 = linspace(2.1,3.2, 10);
    y2 = linspace(-3.3, 3.2, 10);
    plot(y1, 7.5 * ones(10, 1), 'b.', 'LineWidth', 3, 'MarkerSize', 7)
    plot(y1, -12 * ones(10, 1), 'b.', 'LineWidth', 3, 'MarkerSize', 7)
end
% Colors:
%colormap turbo;
%c                   = colorbar('eastoutside'); %, 'Direction','reverse');
%c.Label.String      = '$\left| \Psi  \right|^2$';
%c.Label.Interpreter = 'latex';
%c.Label.FontSize    = FontSize + 5;
% Additional Labels:
%text( -8.7, 65, 0, '(a)', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize);
if Version == 1
    text( 0, 0.45, 0, '$\epsilon > 0$','HorizontalAlignment', 'center', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize + 15);
elseif Version == 2
    text( 3, -6, 0, '$\epsilon > 0$', 'HorizontalAlignment', 'center', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize + 5);
end
% PaperSize:
set(gcf,'paperunits','in');
set(gcf,'papersize',[Position(3) + 1, Position(4) + 1]);
% Saving:
if Version == 1
    fname   = sprintf('Fig_WaveFunction_Potential_NegEps_v1.pdf');
elseif Version == 2
    fname   = sprintf('Fig_WaveFunction_Potential_NegEps_v2.pdf');
end
hfig    = gcf;
print(hfig,'-bestfit','-dpdf', '-r960', fname);

%% Figure 4: Polarization 3D plot theoretical
%
clear all
close all
clc

FontSize = 20;
Position = [1 1 7 6];

Figure = figure('color', 'white', 'units', 'inches', 'Position', Position);

hold on
% Fontsize:
set(gca, 'FontSize', FontSize)

% Axis limits:
xlim([5 9])
%xticks([-7.5 -5 -2.5 0 2.5 5 7.5])
ylim([-0.1 0.1])
yticks([-0.1  0 0.1])
%zlim([])
%axis square
box on
% Axis Labels:
xlabel('$\alpha$', 'Interpreter', 'latex', 'FontSize', FontSize)
ylabel('$\epsilon$', 'FontWeight', 'Bold', 'Interpreter', 'latex', 'FontWeight', 'Bold', 'FontSize', FontSize + 5)
zlabel('$P(\alpha, \epsilon)$', 'Interpreter', 'latex', 'FontSize', FontSize)
% Data:
Polarization = load('Polarization3Particle.mat');
Polarization = Polarization.Polarization;
% Plotting:
for alphaInd = 1:26
    plot3(Polarization.Alpha(alphaInd) * ones(201, 1), Polarization.Kappa, Polarization.Polarization(alphaInd, :, :), 'Color', [alphaInd/26 0.5 (1 - alphaInd/26)], 'LineWidth', 2)
end
view([120 30])
grid
% Colors:
%colormap turbo;
%c                   = colorbar('eastoutside'); %, 'Direction','reverse');
%c.Label.String      = '$\left| \Psi  \right|^2$';
%c.Label.Interpreter = 'latex';
%c.Label.FontSize    = FontSize + 5;
% Additional Labels:
text( 8, -0.172, 3, '(a)', 'HorizontalAlignment', 'center', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize);
%text( 0, 45, 0, '$\epsilon < 0$','HorizontalAlignment', 'center', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize);
% PaperSize:
set(gcf,'paperunits','in');
set(gcf,'papersize',[Position(3) + 1, Position(4) + 1]);
% Saving:
fname   = sprintf('Fig_Polarization_3D_theor.pdf');
hfig    = gcf;
print(hfig,'-bestfit','-dpdf', '-r960', fname);

%% Figure 5: Spectral gaps for 1,3,5 and 7 particles EXPERIMENTAL
%
clear all
close all
clc

FontSize = 20;
Position = [1 1 7 3];

Figure = figure('color', 'white', 'units', 'inches', 'Position', Position);

hold on
% Fontsize:
set(gca, 'FontSize', FontSize)

% Axis limits:
%xlim([0 1100])
%xticks([-200 200 600 1000])
ylim([10^-3 1])
%yticks([0 20 40 60])
%zlim([])
%axis square
%set(gca, 'Yscale', 'log')
%set(gca, 'Xscale', 'log')
%grid
box on
scalingVec(1, :) = [7.182 0.875 0.917 0.482]; % pascu's 1st
scalingVec(2, :) = [13.05 1.327 0.682 0.51];    % Pascu's 2nd
scalingVec(3, :) = [12.8 9.2 9.067 8.02];   % Original from Shahal
scalingVec(4, :) = [9 4 3.5 3.5];   % Eyeballed

Exponent    = 1; %1/2.1; % 1.9 - 2.2 majdnem ugyanúgy szemre jót kapunk.
Factor      = 1;

%xline([114 368 711 1031])
%yline([7.182 0.875 0.917 0.482])
%xline([104 360 720 1020])
%yline([13.05 1.327 0.682 0.51])
%xline([104 280 560 850])
%yline([13.05 9.2 9.067 8.2])
yline(1)
% Axis Labels:
xlabel('$V_k$', 'Interpreter', 'latex', 'FontSize', FontSize)
ylabel('$\Pi^{-1} / \Pi^{-1}_0$', 'FontWeight', 'Bold', 'Interpreter', 'latex', 'FontWeight', 'Bold', 'FontSize', FontSize)
%zlabel('$P(\alpha, \epsilon)$', 'Interpreter', 'latex', 'FontSize', FontSize)
% Data:
Exp1 = load('P1.tsv');
Exp3 = load('P3.tsv');
Exp5 = load('P5.tsv');
Exp7 = load('P7.tsv');
% Plotting:
scaleInd = 3;
%plot(Exp1(:, 1), Exp1(:, 2)/max(Exp1(:, 2)), 'x', 'MarkerSize', 8, 'LineWidth', 1.5)
plot(Factor * Exp1(:, 1).^Exponent, Exp1(:, 2)/scalingVec(scaleInd, 1), 'x', 'MarkerSize', 8, 'LineWidth', 1.5)
%plot(Exp3(:, 1), Exp3(:, 2)/max(Exp3(:, 2)), 'x', 'MarkerSize', 8, 'LineWidth', 1.5)
plot(Factor * Exp3(:, 1).^Exponent, Exp3(:, 2)/scalingVec(scaleInd, 2), 'x', 'MarkerSize', 8, 'LineWidth', 1.5)
%plot(Exp5(:, 1), Exp5(:, 2)/max(Exp5(:, 2)), 'x', 'MarkerSize', 8, 'LineWidth', 1.5)
plot(Factor * Exp5(:, 1).^Exponent, Exp5(:, 2)/scalingVec(scaleInd, 3), 'x', 'MarkerSize', 8, 'LineWidth', 1.5)
%plot(Exp7(:, 1), Exp7(:, 2)/max(Exp7(:, 2)), 'x', 'MarkerSize', 8, 'LineWidth', 1.5)
plot(Factor * Exp7(:, 1).^Exponent, Exp7(:, 2)/scalingVec(scaleInd, 4), 'x', 'MarkerSize', 8, 'LineWidth', 1.5)
% Colors:
%colormap turbo;
%c                   = colorbar('eastoutside'); %, 'Direction','reverse');
%c.Label.String      = '$\left| \Psi  \right|^2$';
%c.Label.Interpreter = 'latex';
%c.Label.FontSize    = FontSize + 5;
% Additional Labels:
%text( 40, 1.7, 0, '(a)', 'HorizontalAlignment', 'center', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize);
%text( 80, 0.5*10^-1, 0, '$1e^-$', 'HorizontalAlignment', 'center', 'FontWeight', 'Bold', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize);
%text( 300, 0.5*10^-1, 0, '$3e^-$', 'HorizontalAlignment', 'center', 'FontWeight', 'Bold', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize);
%text( 600, 0.5*10^-1, 0, '$5e^-$', 'HorizontalAlignment', 'center', 'FontWeight', 'Bold', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize);
%text( 900, 0.5*10^-1, 0, '$7e^-$', 'HorizontalAlignment', 'center', 'FontWeight', 'Bold', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize);

qw{1} = plot(nan, 'kx', 'MarkerSize', 8, 'LineWidth', 1.5)
legend([qw{:}], {'Exp.'}, 'location', 'best')
% PaperSize:
set(gcf,'paperunits','in');
set(gcf,'papersize',[Position(3) + 1, Position(4) + 1]);
% Saving:
fname   = sprintf('Fig_spectral_gap_exp_log_lin.pdf');
hfig    = gcf;
print(hfig,'-bestfit','-dpdf', '-r960', fname);
%% Figure 6: Spectral gaps for 1,3,5 and 7 particles THEORETICAL
%
clear all
%close all
clc

FontSize = 20;
Position = [5 5 7 3];

Figure = figure('color', 'white', 'units', 'inches', 'Position', Position);

hold on
% Fontsize:
set(gca, 'FontSize', FontSize)

% Axis limits:
xlim([0 16])
xticks([0 3 6 9 12 15])
ylim([10^-3 1])
yticks([10^-2 10^-1 1])
%zlim([])
%axis square
set(gca, 'Yscale', 'log')
%set(gca, 'Xscale', 'log')
%grid
box on
% Axis Labels:
xlabel('$\alpha$', 'Interpreter', 'latex', 'FontSize', FontSize)
ylabel('$\Delta / \Delta_0$', 'FontWeight', 'Bold', 'Interpreter', 'latex', 'FontWeight', 'Bold', 'FontSize', FontSize)
%zlabel('$P(\alpha, \epsilon)$', 'Interpreter', 'latex', 'FontSize', FontSize)
% Data:
DMRG1   = load('DMRG_Ne_1_gauss.dat');
DMRG3   = load('DMRG3.mat'); DMRG3 = DMRG3.DMRG3
DMRG5   = load('DMRG5.mat'); DMRG5 = DMRG5.DMRG5
DMRG7   = load('DMRG7.mat'); DMRG7 = DMRG7.DMRG7

ED1 = load('EDSplitting_1_particle_Nx_5000.mat'); ED1 = ED1.data;
ED3 = load('EDSplitting_3_particles_restricted_Nx_90_beta_1e-05.mat'); ED3 = ED3.data;
ED5 = load('ED5P.mat'); ED5 = ED5.data;
ED7 = load('ED7P.mat'); ED7 = ED7.ED7;

IT1 = load('Standard1particleSplitting.mat'); IT1 = IT1.OneParticleInstanton;
IT3 = load('Standar3particledSplitting.mat'); IT3 = IT3.SPLITTINGS;
IT5 = load('Standard5particleSplitting.mat'); IT5 = IT5.SPLITTINGS;
IT7e = load('Standard7particleSplitting.mat');
IT7 = IT7e.SPLITTINGS; %IT7e.Instanton7;

% Plotting:
xshift = 0;
xscale = 1;
yscale = 1;
xline(6.8)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot((ED1(:, 1) + xshift) * xscale, ED1(:, 2) * yscale, 'rs', 'MarkerSize', 10, 'LineWidth', 1.5)
plot((-DMRG3(10:end-9, 3) + xshift) * xscale, abs(DMRG3(10:end-9, 6)) * yscale, 'rs', 'MarkerSize', 10, 'LineWidth', 1.5)
plot((-DMRG5(1:end-17, 3) + xshift) * xscale + 0.05, DMRG5(1:end-17, 6) * yscale, 'rs', 'MarkerSize', 10, 'LineWidth', 1.5)
plot((DMRG7(1, 1:end-2) + xshift) * xscale, DMRG7(2, 1:end-2) * yscale, 'rs', 'MarkerSize', 10, 'LineWidth', 1.5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot((ED1(:, 1) + xshift) * xscale, ED1(:, 2) * yscale, 'b.-', 'MarkerSize', 15)
plot((ED3(:, 1) + xshift) * xscale, ED3(:, 2) * yscale, 'b.-', 'MarkerSize', 15)
plot((ED5(:, 1) + xshift) * xscale, ED5(:, 2) * yscale, 'b.-', 'MarkerSize', 15)
plot((ED7(1:end-2, 1) - 0.05 + xshift) * xscale, ED7(1:end-2, 2) * yscale, 'b.-', 'MarkerSize', 15)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot((IT1(:, 1) + xshift) * xscale, IT1(:, 2) * yscale, 'ko--', 'MarkerFaceColor', 'k')
plot((IT3(:, 1) + xshift) * xscale, IT3(:, 3) * yscale, 'ko--', 'MarkerFaceColor', 'k')
plot((IT5(:, 1) + xshift) * xscale, IT5(:, 3) * yscale, 'ko--', 'MarkerFaceColor', 'k')
plot((IT7(:, 1) + xshift) * xscale, IT7(:, 2) * yscale, 'ko--', 'MarkerFaceColor', 'k')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Colors:
%colormap turbo;
%c                   = colorbar('eastoutside'); %, 'Direction','reverse');
%c.Label.String      = '$\left| \Psi  \right|^2$';
%c.Label.Interpreter = 'latex';
%c.Label.FontSize    = FontSize + 5;
% Additional Labels:
text( 0.6, 2, 0, '(b)', 'HorizontalAlignment', 'center', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize);
%text( 0, 45, 0, '$\epsilon < 0$','HorizontalAlignment', 'center', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize);
% Legend:


qw{1} = plot(nan, 'ko--', 'MarkerFaceColor', 'k');
qw{2} = plot(nan, 'b.-', 'MarkerSize', 15);
qw{3} = plot(nan, 'rs', 'MarkerSize', 10, 'LineWidth', 1.5)
%qw{4} = plot(nan, 'k-d'); % You can add an extra element too
legend([qw{:}], {'IT','ED','DMRG'}, 'Position', [0.4 3.2*10^-1 1 1])
% PaperSize:
set(gcf,'paperunits','in');
set(gcf,'papersize',[Position(3) + 1, Position(4) + 1]);
% Saving:
fname   = sprintf('Fig_spectral_gap_theor.pdf');
hfig    = gcf;
print(hfig,'-bestfit','-dpdf', '-r960', fname);
%% Figure ??: Boxes
clear all
close all
clc

FontSize = 20;
Position = [5 5 7 3];

Figure = figure('color', 'white', 'units', 'inches', 'Position', Position);

hold on
% Fontsize:
set(gca, 'FontSize', FontSize)

% Axis limits:
xlim([0 16])
xticks([0 3 6 9 12 15])
ylim([10^-3 1])
yticks([10^-2 10^-1 1])
%zlim([])
%axis square
set(gca, 'Yscale', 'log')
%set(gca, 'Xscale', 'log')
%grid
box on
% Axis Labels:
xlabel('$\alpha$', 'Interpreter', 'latex', 'FontSize', FontSize)
ylabel('$\Delta / \Delta_0$', 'FontWeight', 'Bold', 'Interpreter', 'latex', 'FontWeight', 'Bold', 'FontSize', FontSize)
%zlabel('$P(\alpha, \epsilon)$', 'Interpreter', 'latex', 'FontSize', FontSize)
% Data:
xshift = 0;
xscale = 1;
yscale = 1/1;
DMRG1   = load('DMRG_Ne_1_gauss.dat');
DMRG3   = load('DMRG3.mat'); DMRG3 = DMRG3.DMRG3
DMRG5   = load('DMRG5.mat'); DMRG5 = DMRG5.DMRG5
DMRG7   = load('DMRG7.mat'); DMRG7 = DMRG7.DMRG7

ED1 = load('EDSplitting_1_particle_Nx_5000.mat'); ED1 = ED1.data;
ED3 = load('EDSplitting_3_particles_restricted_Nx_90_beta_1e-05.mat'); ED3 = ED3.data;
ED5 = load('ED5P.mat'); ED5 = ED5.data;
ED7 = load('ED7P.mat'); ED7 = ED7.ED7;

IT1 = load('Standard1particleSplitting.mat'); IT1 = IT1.OneParticleInstanton;
IT3 = load('Standar3particledSplitting.mat'); IT3 = IT3.SPLITTINGS;
IT5 = load('Standard5particleSplitting.mat'); IT5 = IT5.SPLITTINGS;
IT7e = load('Standard7particleSplitting.mat');
IT7 = IT7e.SPLITTINGS; %IT7e.Instanton7;

Exp1 = load('P1.tsv');
Exp3 = load('P3.tsv');
Exp5 = load('P5.tsv');
Exp7 = load('P7.tsv');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot((ED1(:, 1) + xshift) * xscale, ED1(:, 2) * yscale, 'rs', 'MarkerSize', 10, 'LineWidth', 1.5)
plot((-DMRG3(10:end-9, 3) + xshift) * xscale, abs(DMRG3(10:end-9, 6)) * yscale, 'rs', 'MarkerSize', 10, 'LineWidth', 1.5)
plot((-DMRG5(1:end-17, 3) + xshift) * xscale + 0.05, DMRG5(1:end-17, 6) * yscale, 'rs', 'MarkerSize', 10, 'LineWidth', 1.5)
plot((DMRG7(1, 1:end-2) + xshift) * xscale, DMRG7(2, 1:end-2) * yscale, 'rs', 'MarkerSize', 10, 'LineWidth', 1.5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot((ED1(:, 1) + xshift) * xscale, ED1(:, 2) * yscale, 'b.-', 'MarkerSize', 15)
plot((ED3(:, 1) + xshift) * xscale, ED3(:, 2) * yscale, 'b.-', 'MarkerSize', 15)
plot((ED5(:, 1) + xshift) * xscale, ED5(:, 2) * yscale, 'b.-', 'MarkerSize', 15)
plot((ED7(1:end-2, 1) - 0.05 + xshift) * xscale, ED7(1:end-2, 2) * yscale, 'b.-', 'MarkerSize', 15)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot((IT1(:, 1) + xshift) * xscale, IT1(:, 2) * yscale, 'ko--', 'MarkerFaceColor', 'k')
plot((IT3(:, 1) + xshift) * xscale, IT3(:, 3) * yscale, 'ko--', 'MarkerFaceColor', 'k')
plot((IT5(:, 1) + xshift) * xscale, IT5(:, 3) * yscale, 'ko--', 'MarkerFaceColor', 'k')
plot((IT7(:, 1) + xshift) * xscale, IT7(:, 2) * yscale, 'ko--', 'MarkerFaceColor', 'k')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Figure 7: Spectral gap for 3 BOTH
%
clear all
close all
clc

FontSize = 20;
Position = [1 1 7 6];

Figure = figure('color', 'white', 'units', 'inches', 'Position', Position);

hold on
% Fontsize:
set(gca, 'FontSize', FontSize)

% Axis limits:
xlim([4.5 8.5])
%xticks([0 3 6 9 12 15])
ylim([10^-2 1])
%yticks([10^-2 10^-1 1])
%zlim([])
%axis square
set(gca, 'Yscale', 'log')
%grid
box on
% Axis Labels:
xlabel('$\alpha$', 'Interpreter', 'latex', 'FontSize', FontSize + 5)
ylabel('$\Delta$', 'FontWeight', 'Bold', 'Interpreter', 'latex', 'FontWeight', 'Bold', 'FontSize', FontSize + 5)
%zlabel('$P(\alpha, \epsilon)$', 'Interpreter', 'latex', 'FontSize', FontSize)
% Data:
Exp3 = load('P3.tsv');

DMRG3   = load('DMRG3.mat'); DMRG3 = DMRG3.DMRG3;

ED3 = load('EDSplitting_3_particles_restricted_Nx_90_beta_1e-05.mat'); ED3 = ED3.data;

IT3 = load('Standar3particledSplitting.mat'); IT3 = IT3.SPLITTINGS;

% Plotting:
text( 4.75, 0.75, 0, '(c)', 'HorizontalAlignment', 'center', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize);
axis square
xshift = -20;
xscale = 1/52;
yscale = 1/10;
plot((Exp3(:, 1) - xshift) * xscale, Exp3(:, 2) * yscale, 'x', 'Color', [0.8500 0.3250 0.0980], 'MarkerSize', 10, 'LineWidth', 1.5)
xshift = 0;
xscale = 1;
yscale = 1;
plot((-DMRG3(10:end-9, 3) + xshift) * xscale, abs(DMRG3(10:end-9, 6)) * yscale, 'rs', 'MarkerSize', 10, 'LineWidth', 1.5)
plot((ED3(:, 1) + xshift) * xscale, ED3(:, 2) * yscale, 'b.-', 'MarkerSize', 15)
plot((IT3(:, 1) + xshift) * xscale, IT3(:, 3) * yscale, 'ko--', 'MarkerFaceColor', 'k')
hold off
axes('Position',[.29 .26 .34 .34])
box on
hold on
scatter([2.85 7.4 10.8 13.7], [114 368 711 1031], 100, 'r', 'Filled')
scatter([7.4], [368], 150, 'bd', 'Filled')
x = linspace(0, 15, 1000);
plot(x, (x -2.5)*90, 'k--', 'LineWidth', 4)
xlim([1.5 14.5])
ylim([0 1100])
ax = gca;
ax.FontSize = 20;
xlabel('\alpha_Q', 'FontSize', 20, 'Position', [15 -30])
ylabel('V_Q [V]', 'FontSize', 20)
xticks([3 6 9 12])
yticks([0 200 400 600 800 1000])
yticklabels({'0.0','0.2','0.4','0.6','0.8','1'})
axis square
hold off
% Colors:
%colormap turbo;
%c                   = colorbar('eastoutside'); %, 'Direction','reverse');
%c.Label.String      = '$\left| \Psi  \right|^2$';
%c.Label.Interpreter = 'latex';
%c.Label.FontSize    = FontSize + 5;
% Additional Labels:

%text( 0, 45, 0, '$\epsilon < 0$','HorizontalAlignment', 'center', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize);
% Legend:

% PaperSize:
set(gcf,'paperunits','in');
set(gcf,'papersize',[Position(3) + 1, Position(4) + 1]);
% Saving:
fname   = sprintf('Fig_spectral_gap_comparison_N_3.pdf');
hfig    = gcf;
print(hfig,'-bestfit','-dpdf', '-r960', fname);
%% Figure 7: Spectral gap for 1 BOTH
%
clear all
close all
clc

FontSize = 20;
Position = [1 1 7 6];

Figure = figure('color', 'white', 'units', 'inches', 'Position', Position);

hold on
% Fontsize:
set(gca, 'FontSize', FontSize)

% Axis limits:
xlim([0 5])
%xticks([0 3 6 9 12 15])
ylim([10^-2 1])
%yticks([10^-2 10^-1 1])
%zlim([])
%axis square
set(gca, 'Yscale', 'log')
%grid
box on
% Axis Labels:
xlabel('$\alpha$', 'Interpreter', 'latex', 'FontSize', FontSize + 5)
ylabel('$\Delta$', 'FontWeight', 'Bold', 'Interpreter', 'latex', 'FontWeight', 'Bold', 'FontSize', FontSize + 5)
%zlabel('$P(\alpha, \epsilon)$', 'Interpreter', 'latex', 'FontSize', FontSize)
% Data:
Exp1 = load('P1.tsv');



ED1 = load('EDSplitting_1_particle_Nx_5000.mat'); ED1 = ED1.data;
DMRG1   = ED1;
IT1 = load('Standard1particleSplitting.mat'); IT1 = IT1.OneParticleInstanton;

% Plotting:
text( 4.75, 0.75, 0, '(c)', 'HorizontalAlignment', 'center', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize);
axis square
xshift = 65;
xscale = 1/20;
yscale = 1/30;
plot((Exp1(1:end-4, 1) - xshift) * xscale, Exp1(1:end-4, 2) * yscale, 'x', 'Color', [0.8500 0.3250 0.0980], 'MarkerSize', 10, 'LineWidth', 1.5)
xshift = 0;
xscale = 1;
yscale = 1;
plot((ED1(:, 1) + xshift) * xscale, ED1(:, 2) * yscale, 'rs', 'MarkerSize', 10, 'LineWidth', 1.5)
plot((ED1(:, 1) + xshift) * xscale, ED1(:, 2) * yscale, 'b.-', 'MarkerSize', 15)
plot((IT1(:, 1) + xshift) * xscale, IT1(:, 2) * yscale, 'ko--', 'MarkerFaceColor', 'k')
hold off
%axes('Position',[.29 .26 .34 .34])
%box on
%hold on
%scatter([2.85 7.4 10.8 13.7], [114 368 711 1031], 100, 'r', 'Filled')
%scatter([7.4], [368], 150, 'bd', 'Filled')
x = linspace(0, 15, 1000);
%plot(x, (x -2.5)*90, 'k--', 'LineWidth', 4)
%xlim([1.5 14.5])
%ylim([0 1100])
%ax = gca;
%ax.FontSize = 20;
%xlabel('\alpha_Q', 'FontSize', 20, 'Position', [15 -30])
%ylabel('V_Q [V]', 'FontSize', 20)
%xticks([3 6 9 12])
%yticks([0 200 400 600 800 1000])
%yticklabels({'0.0','0.2','0.4','0.6','0.8','1'})
axis square
hold off
% Colors:
%colormap turbo;
%c                   = colorbar('eastoutside'); %, 'Direction','reverse');
%c.Label.String      = '$\left| \Psi  \right|^2$';
%c.Label.Interpreter = 'latex';
%c.Label.FontSize    = FontSize + 5;
% Additional Labels:

%text( 0, 45, 0, '$\epsilon < 0$','HorizontalAlignment', 'center', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize);
% Legend:

% PaperSize:
set(gcf,'paperunits','in');
set(gcf,'papersize',[Position(3) + 1, Position(4) + 1]);
% Saving:
fname   = sprintf('Fig_spectral_gap_comparison_N_1.pdf');
hfig    = gcf;
print(hfig,'-bestfit','-dpdf', '-r960', fname);

%% Figure 7: Spectral gap for 5 BOTH
%
clear all
close all
clc

FontSize = 20;
Position = [1 1 7 6];

Figure = figure('color', 'white', 'units', 'inches', 'Position', Position);

hold on
% Fontsize:
set(gca, 'FontSize', FontSize)

% Axis limits:
xlim([9 12])
%xticks([0 3 6 9 12 15])
ylim([10^-2 1])
%yticks([10^-2 10^-1 1])
%zlim([])
%axis square
set(gca, 'Yscale', 'log')
%grid
box on
% Axis Labels:
xlabel('$\alpha$', 'Interpreter', 'latex', 'FontSize', FontSize + 5)
ylabel('$\Delta$', 'FontWeight', 'Bold', 'Interpreter', 'latex', 'FontWeight', 'Bold', 'FontSize', FontSize + 5)
%zlabel('$P(\alpha, \epsilon)$', 'Interpreter', 'latex', 'FontSize', FontSize)
% Data:
Exp5 = load('P5.tsv');



ED5 = load('ED5P.mat'); ED5 = ED5.data;
DMRG5   = load('DMRG5.mat'); DMRG5 = DMRG5.DMRG5
IT5 = load('Standard5particleSplitting.mat'); IT5 = IT5.SPLITTINGS;

% Plotting:
text( 4.75, 0.75, 0, '(c)', 'HorizontalAlignment', 'center', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize);
axis square
xshift = -490;
xscale = 1/110;
yscale = 1/10;
plot((Exp5(1:end-4, 1) - xshift) * xscale, Exp5(1:end-4, 2) * yscale, 'x', 'Color', [0.8500 0.3250 0.0980], 'MarkerSize', 10, 'LineWidth', 1.5)
xshift = 0;
xscale = 1;
yscale = 1;
plot((ED5(:, 1) + xshift) * xscale, ED5(:, 2) * yscale, 'b.-', 'MarkerSize', 15)

plot((-DMRG5(1:end-17, 3) + xshift) * xscale + 0.05, DMRG5(1:end-17, 6) * yscale, 'rs', 'MarkerSize', 10, 'LineWidth', 1.5)
plot((IT5(:, 1) + xshift) * xscale, IT5(:, 3) * yscale, 'ko--', 'MarkerFaceColor', 'k')
hold off
%axes('Position',[.29 .26 .34 .34])
%box on
%hold on
%scatter([2.85 7.4 10.8 13.7], [114 368 711 1031], 100, 'r', 'Filled')
%scatter([7.4], [368], 150, 'bd', 'Filled')
x = linspace(0, 15, 1000);
%plot(x, (x -2.5)*90, 'k--', 'LineWidth', 4)
%xlim([1.5 14.5])
%ylim([0 1100])
%ax = gca;
%ax.FontSize = 20;
%xlabel('\alpha_Q', 'FontSize', 20, 'Position', [15 -30])
%ylabel('V_Q [V]', 'FontSize', 20)
%xticks([3 6 9 12])
%yticks([0 200 400 600 800 1000])
%yticklabels({'0.0','0.2','0.4','0.6','0.8','1'})
axis square
hold off
% Colors:
%colormap turbo;
%c                   = colorbar('eastoutside'); %, 'Direction','reverse');
%c.Label.String      = '$\left| \Psi  \right|^2$';
%c.Label.Interpreter = 'latex';
%c.Label.FontSize    = FontSize + 5;
% Additional Labels:

%text( 0, 45, 0, '$\epsilon < 0$','HorizontalAlignment', 'center', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize);
% Legend:

% PaperSize:
set(gcf,'paperunits','in');
set(gcf,'papersize',[Position(3) + 1, Position(4) + 1]);
% Saving:
fname   = sprintf('Fig_spectral_gap_comparison_N_5.pdf');
hfig    = gcf;
print(hfig,'-bestfit','-dpdf', '-r960', fname);
%% Figure 8: Perp factor
%
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
xlim([1 16])
xticks([0 3 6 9 12 15])
ylim([0 40])
yticks([0 20 40])
%zlim([])
%axis square
%set(gca, 'Yscale', 'log')
%grid
box on
% Axis Labels:
xlabel('$\alpha$', 'Interpreter', 'latex', 'FontSize', FontSize)
ylabel('$R_0(\alpha, N)$', 'FontWeight', 'Bold', 'Interpreter', 'latex', 'FontWeight', 'Bold', 'FontSize', FontSize)
%zlabel('$P(\alpha, \epsilon)$', 'Interpreter', 'latex', 'FontSize', FontSize)
% Data:
IT1 = load('Standard1particleSplitting.mat'); IT1 = IT1.OneParticleInstanton;
IT3 = load('Standar3particledSplitting.mat'); IT3 = IT3.SPLITTINGS;
IT5 = load('Standard5particleSplitting.mat'); IT5 = IT5.SPLITTINGS;
IT7P= load('IT7Perp.mat'); IT7P = IT7P.IT7;
% Plotting:
plot(0:0.5:8, ones(1, length(0:0.5:8)), '-', 'LineWidth', 3)
plot(IT3(:, 1), IT3(:, 3) ./ IT3(:, 2), '-', 'LineWidth', 3)
plot(IT5(:, 1), IT5(:, 3) ./ IT5(:, 2), '-', 'LineWidth', 3)
plot(IT7P(1,:), IT7P(2, :), '-', 'LineWidth', 3)

scatter([2.9 7.4 10.8 13.7], [1 3.1 7.58 23], 100, 'rd', 'Filled')
% Colors:
%colormap turbo;
%c                   = colorbar('eastoutside'); %, 'Direction','reverse');
%c.Label.String      = '$\left| \Psi  \right|^2$';
%c.Label.Interpreter = 'latex';
%c.Label.FontSize    = FontSize + 5;
% Additional Labels:
%text( 2, 36, 3, '(c)', 'HorizontalAlignment', 'center', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize);
%text( 0, 45, 0, '$\epsilon < 0$','HorizontalAlignment', 'center', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize);
text( 2.9, 5, 0, '$1e^-$', 'HorizontalAlignment', 'center', 'FontWeight', 'Bold', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize);
text( 7.4, 8, 0, '$3e^-$', 'HorizontalAlignment', 'center', 'FontWeight', 'Bold', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize);
text( 10.8, 13, 0, '$5e^-$', 'HorizontalAlignment', 'center', 'FontWeight', 'Bold', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize);
text( 13.2, 25, 0, '$7e^-$', 'HorizontalAlignment', 'center', 'FontWeight', 'Bold', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize);

% PaperSize:
set(gcf,'paperunits','in');
set(gcf,'papersize',[Position(3) + 0.2, Position(4) + 0.2]);
% Saving:
fname   = sprintf('Fig_perpfactors.pdf');
hfig    = gcf;
print(hfig,'-bestfit','-dpdf', '-r960', fname);

