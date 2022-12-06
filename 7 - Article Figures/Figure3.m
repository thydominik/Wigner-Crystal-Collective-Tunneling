clc
clear all

% Experimental data, unscaled

Exp1 = load('P1.tsv');
Exp3 = load('P3.tsv');
Exp5 = load('P5.tsv');
Exp7 = load('P7.tsv');

figure(41)
clf(figure(41))
hold on
scale = 1/200;
plot(Exp1(:, 1) * scale, Exp1(:, 2))
plot(Exp3(:, 1) * scale, Exp3(:, 2))
plot(Exp5(:, 1) * scale, Exp5(:, 2))
plot(Exp7(:, 1) * scale, Exp7(:, 2))
hold off
% DMRG:

DMRG1   = load('DMRG_Ne_1_gauss.dat');
DMRG3   = load('DMRG3.mat'); DMRG3 = DMRG3.DMRG3
DMRG5   = load('DMRG5.mat'); DMRG5 = DMRG5.DMRG5
DMRG7   = load('DMRG7.mat'); DMRG7 = DMRG7.DMRG7


figure(42)
clf(figure(42))
hold on
plot(-DMRG1(:, 3), abs(DMRG1(:, 4) - DMRG1(:, 5)))
plot(-DMRG3(:, 3), abs(DMRG3(:, 4) - DMRG3(:, 5)))
plot(-DMRG5(:, 3), abs(DMRG5(:, 4) - DMRG5(:, 5)))
plot(DMRG7(1, :), DMRG7(2, :))
hold off
% Exact diagonalization:

ED1 = load('EDSplitting_1_particle_Nx_5000.mat'); ED1 = ED1.data;
ED3 = load('EDSplitting_3_particles_restricted_Nx_90_beta_1e-05.mat'); ED3 = ED3.data;
ED5 = load('ED5P.mat'); ED5 = ED5.data;
ED7 = load('ED7P.mat'); ED7 = ED7.ED7;

figure(43)
clf(figure(43))
hold on
plot(ED1(:, 1), ED1(:, 2))
plot(ED3(:, 1), ED3(:, 2))
plot(ED5(:, 1), ED5(:, 2))
plot(ED7(:, 1), ED7(:, 2))
hold off

% Instanton:

IT1 = load('Standard1particleSplitting.mat'); IT1 = IT1.OneParticleInstanton;
IT3 = load('Standar3particledSplitting.mat'); IT3 = IT3.SPLITTINGS;
IT5 = load('Standard5particleSplitting.mat'); IT5 = IT5.SPLITTINGS;
IT7e = load('Standard7particleSplitting.mat');
IT7 = IT7e.SPLITTINGS; %IT7e.Instanton7;

figure(44)
clf(figure(44))
hold on
plot(IT1(:, 1), IT1(:, 2), '.-')
plot(IT3(:, 1), IT3(:, 3), '.-')
plot(IT5(:, 1), IT5(:, 3), '.-')
plot(IT7(:, 1), IT7(:, 2), '.-')
hold off
%
figure(4)
clf(figure(4))
figure('Renderer', 'painters', 'Position', [10 10 900 600])
t = tiledlayout(2, 1);

ax1 = nexttile;
ylabel('\Delta_{Exp}', 'FontSize', 20)
set(ax1,'xticklabel',[])
grid
hold on
ylim([2*10^-1 0.3*10^2])
xlim([0 18])
box
xshift = 49/56;
xscale = 1/56;
yscale = 1
plot((Exp1(:, 1) - xshift*90) * xscale * 3, Exp1(:, 2) * yscale/2, 'x', 'MarkerSize', 8, 'LineWidth', 1.5)
plot((Exp3(:, 1) - xshift) * xscale, Exp3(:, 2) * yscale, 'x', 'MarkerSize', 8, 'LineWidth', 1.5)
plot((Exp5(:, 1) - xshift) * xscale, Exp5(:, 2) * yscale, 'x', 'MarkerSize', 8, 'LineWidth', 1.5)
plot((Exp7(:, 1) - xshift) * xscale, Exp7(:, 2) * yscale, 'x', 'MarkerSize', 8, 'LineWidth', 1.5)
set(gca,'Yscale', 'log')
qw{1} = plot(nan, 'kx', 'MarkerSize', 8, 'LineWidth', 1.5)
legend([qw{:}], {'experimental data'}, 'location', 'best')
hold off

ax2 = nexttile;
ylabel('\Delta_{Th.}', 'FontSize', 20)
xlabel('\alpha', 'FontSize', 20)
hold on
box

xshift = 0;
xscale = 1;
yscale = 1;
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
ylim([10^-3 1.5])
set(gca,'Yscale', 'log')
grid
qw{1} = plot(nan, 'ko--', 'MarkerFaceColor', 'k');
qw{2} = plot(nan, 'b.-', 'MarkerSize', 15);
qw{3} = plot(nan, 'rs', 'MarkerSize', 10, 'LineWidth', 1.5)
%qw{4} = plot(nan, 'k-d'); % You can add an extra element too
legend([qw{:}], {'instanton theory','egzact diagonalization','DMRG'}, 'location', 'best')
hold off


%xlabel(t, '\alpha', 'FontSize', 20)
%ylabel(t, '\Delta', 'FontSize', 20)
t.TileSpacing = 'none'
ax1.FontSize = 20
ax2.FontSize = 20
linkaxes([ax1, ax2], 'x')
%%

figure(5)
clf(figure(5))
figure('Renderer', 'painters', 'Position', [10 10 1000 1000])
hold on
ax = gca;
ax.FontSize = 20;
axis square
box
grid
xlabel('\alpha', 'FontSize', 25)
ylabel('\Delta_{Th.} & \Delta_{Scaled Exp.}', 'FontSize', 25)
ylim([10^-2 1.5])
set(gca, 'Yscale', 'log')
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
axes('Position',[.28 .22 .34 .34])
box on
hold on
scatter([2.9 7.4 10.8 13.7], [114 368 711 1031], 100, 'r', 'Filled')
scatter([7.4], [368], 150, 'bd', 'Filled')
x = linspace(0, 15, 1000);
plot(x, (x -2.5)*90, 'k--', 'LineWidth', 4)
xlim([1.5 14.5])
ylim([0 1100])
ax = gca;
ax.FontSize = 20;
xlabel('\alpha_Q', 'FontSize',20)
ylabel('V_Q [V]', 'FontSize', 20)
xticks([3 6 9 12])
yticks([0 200 400 600 800 1000])
yticklabels({'0.0','0.2','0.4','0.6','0.8','1'})
axis square
hold off

%%
figure(6)
clf(figure(6))
hold on
scatter([2.9 7.4 10.8 13.7], [114 368 711 1031], 100, 'r', 'Filled')
x = linspace(0, 15, 1000);
plot(x, (x -2.5)*90, 'k--', 'LineWidth', 4)
xlim([1.5 14.5])
ylim([0 1100])
ax = gca;
ax.FontSize = 50;
xlabel('\alpha_Q', 'FontSize', 40)
ylabel('V_Q [mV]', 'FontSize', 40)
axis square
hold off



