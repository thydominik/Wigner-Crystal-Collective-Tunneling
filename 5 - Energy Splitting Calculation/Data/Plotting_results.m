clc
clear all

D = load('Traj_3p_STDMC9.mat');
D = D.IterData;
figure(16)
clf(figure(16))
hold on
plot(D.time, D.Trajectories, 'LineWidth', 4)
ylabel('\chi_i', 'FontSize', 20)
xlabel('z', 'FontSize', 20)
yline(0)
xline(0)
hold off


%%
clc
clear all

Experimental = 0;
ExpD1 = load('P1.tsv');
ExpD3 = load('P3.tsv');
ExpD5 = load('P5.tsv');
ExpD7 = load('P7.tsv');

Instanton = 1;
D1 = load('OneParticleInstanton.mat');
D3 = load('StandardSplitting.mat'); D3 = D3.SPLITTINGS;
D5 = load('Standard5particleSplitting.mat'); D5 = D5.SPLITTINGS;

ED = 1;
ED1 = load('EDSplitting_1_particle_Nx_5000.mat'); ED1 = ED1.data;
ED3 = load('EDSplitting_3_particles_restricted_Nx_80_beta_1e-05.mat'); ED3 = ED3.data;
ED5 = load('EDSplitting_5_particles_restricted_Nx1_15_Nx2_20_Nx3_40_Nx4_20_Nx5_15_beta_1e-05.mat'); ED5 = ED5.data;

DMRG = 1;
DMRG1 = load('OneParticleDMRG.txt');
DMRG3 = load('Delta_E_DMRG_Norb_8_eta_20.00.mat');
DMRG5 = load('5particleDMRG.txt');

figure(17)
clf(figure(17))
hold on
if Experimental == 1
    %plot((ExpD1(:, 1) - 0),ExpD1(:, 2), 'r*')
    %plot(ExpD3(:, 1), ExpD3(:, 2), 'r*')
    %plot(ExpD5(:, 1), ExpD5(:, 2), 'r*')
    plot(ExpD7(:, 1), ExpD7(:, 2), 'r*')
end

if Instanton == 1
    plot(D1.OneParticleInstanton(1:end, 1), D1.OneParticleInstanton(1:end, 2), 'k.--', 'MarkerSize', 23, 'LineWidth', 2)
    plot(D3(1:end, 1), D3(1:end, 3), 'k.--', 'MarkerSize', 23, 'LineWidth', 2)
    plot(D5(1:end, 1), D5(1:end, 3), 'k.--', 'MarkerSize', 23, 'LineWidth', 2)
end
int = interp1(ED5(:, 1),  ED5(:, 2), linspace(ED5(1,1), ED5(end,1), 100), 'spline');
if ED == 1
    plot(ED1(:, 1), ED1(:, 2), 'color', [0, 0.4470, 0.7410], 'LineWidth', 2)
    plot(ED3(:, 1), ED3(:, 2), 'color', [0, 0.4470, 0.7410], 'LineWidth', 2)
    %plot(ED5(:, 1), ED5(:, 2), 'color', [0, 0.4470, 0.7410], 'LineWidth', 2)
    plot(linspace(ED5(1,1), ED5(end,1), 100), int(:), 'color', [0, 0.4470, 0.7410], 'LineWidth', 2)
end
int2 = interp1(-DMRG5(1:end-2,3) + 0.2,  DMRG5(1:end-2,end), linspace(-DMRG5(1,3) + 0.2, -DMRG5(end-2,3) + 0.2, 50), 'spline');
if DMRG == 1
    plot(-DMRG1(1:end-3, 1), DMRG1(1:end-3, 2), 's', 'LineWidth', 1.5, 'MarkerSize', 7, 'color', [0.8500, 0.3250, 0.0980]	)
    plot(-DMRG3.alpha_list(35:2:end-16), DMRG3.delta_E_DMRG(35:2:end-16), 's', 'LineWidth', 1.5, 'MarkerSize', 7, 'color', [0.8500, 0.3250, 0.0980]	)
    %plot(-DMRG5(1:end-2,3) + 0.2, DMRG5(1:end-2,end), 's', 'LineWidth', 1.5, 'MarkerSize', 7, 'color', [0.6350, 0.0780, 0.1840])
    plot(linspace(-DMRG5(1,3) + 0.2, -DMRG5(end-2,3) + 0.2, 50), int2, 's', 'LineWidth', 1.5, 'MarkerSize', 7, 'color', [0.8500, 0.3250, 0.0980]	)
end



set(gca, 'Yscale', 'log')

xlabel('\alpha', 'FontSize', 20)
ylabel('\Delta', 'FontSize', 20)
xlim([0 14])
ylim([10^-4 1])
hold off