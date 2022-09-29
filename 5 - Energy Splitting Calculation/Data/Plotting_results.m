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
D7 = load('Standard7particleSplitting.mat'); D7 = D7.SPLITTINGS;
D72= load('Standard7particleSplitting2.mat'); D72 = D72.SPLITTINGS;
Instanton7(:, 1) = [D72(3, 1) D7([2 3 5 6 7 8 9], 1).'];
Instanton7(:, 2) = [D72(3, 3) D7([2 3 5 6 7 8 9], 3).'];

ED = 1;
ED1 = load('EDSplitting_1_particle_Nx_5000.mat'); ED1 = ED1.data;
ED3 = load('EDSplitting_3_particles_restricted_Nx_80_beta_1e-05.mat'); ED3 = ED3.data;
ED5 = load('EDSplitting_5_particles_restricted_Nx1_15_Nx2_20_Nx3_40_Nx4_20_Nx5_15_beta_1e-05.mat'); ED5 = ED5.data;
ED7 = load('EDSplitting_7_particles_restricted_Nx1_5_Nx2_10_Nx3_10_Nx4_20_Nx5_10_Nx6_10_Nx7_5_beta_0_01.mat'); ED7 = ED7.data;
ED7B(1, :) = [11 11.666 12.333 13 13.666];
ED7B(2, :) = [2.23 1.67 1.107 0.5034 0.07853427];
temp = load('alpha.mat')
clear ED7
ED7(:, 1) = temp.alpha(1:18);
temp = load('dE.mat');
ED7(:, 2) = temp.dE;





DMRG = 1;
DMRG1 = load('OneParticleDMRG.txt');
DMRG3 = load('Delta_E_DMRG_Norb_8_eta_20.00.mat');
DMRG5 = load('5particleDMRG.txt');
DMRG7 = load('DMRG7Particle.mat');

TEMPD7(:, 1) = [D72(2, 1); D72(3, 1); D7(4:8, 1)];
TEMPD7(:, 2) = [D72(2, 3); D72(3, 3); D7(4:8, 3)];
NInt = 12;
D77(1, :) = interp1(TEMPD7(:, 1), TEMPD7(:, 2), linspace(12,17.5, NInt));
D77(2, :) = linspace(12,17.5, NInt)


figure(17)
clf(figure(17))
hold on
if Experimental == 0
    %plot((ExpD1(:, 1) - 0),ExpD1(:, 2), 'r*')
    %plot(ExpD3(:, 1), ExpD3(:, 2), 'r*')
    %plot(ExpD5(:, 1), ExpD5(:, 2), 'r*')
    %plot(ExpD7(:, 1), ExpD7(:, 2), 'r*')
end

D7Int(1, :) = interp1(Instanton7(:, 1), Instanton7(:, 2), linspace(12,17.5, 10));
D7Int(2, :) = linspace(12,17.5, 10);
NNN = 4;
D7Int(1, NNN) = D7Int(1, NNN) + 0.35 * 10^-1;
D7(1, 9) = D7(1,9) - 2;
D7(3, 9) = D7(3,9) - 0.5;
if Instanton == 1
    %plot(D1.OneParticleInstanton(1:end, 1), D1.OneParticleInstanton(1:end, 2), 'k.--', 'MarkerSize', 23, 'LineWidth', 2)
    plot(D1.OneParticleInstanton(1:end, 1), ones(length(D1.OneParticleInstanton(1:end, 1)),1), 'LineWidth', 2)
    plot(D3(1:end, 1) - 4.45, D3(1:end, 3) ./ D3(1:end, 2), 'LineWidth', 2)
    %plot(D5(1:end, 1), D5(1:end, 3), 'k.--', 'MarkerSize', 23, 'LineWidth', 2)
    plot(D5(1:end, 1) - 7.8, D5(1:end, 3) ./ D5(1:end, 2), 'LineWidth', 2)
    %plot(Instanton7(:, 1), Instanton7(:, 2), 'k.--', 'MarkerSize', 23, 'LineWidth', 2)
    %plot(D7Int(2, :), D7Int(1, :), 'k.--', 'MarkerSize', 23, 'LineWidth', 2)
    plot(D7([1 2 3 4 6 7 8 9], 1) - 10.45, D7([1 2 3 4 6 7 8 9], 9), 'LineWidth', 2)
    %plot(D7([1 2 3 5 6 7 8 9], 1), D7([1 2 3 5 6 7 8 9], 3), 'k.--', 'MarkerSize', 23, 'LineWidth', 2, 'DisplayName', 'Instanton')
    %plot(D72(1:end, 1), D72(1:end, 3), 'k.--', 'MarkerSize', 23, 'LineWidth', 2)
    %plot(TEMPD7(1:end, 1), TEMPD7(1:end, 2), 'k.--', 'MarkerSize', 23, 'LineWidth', 2)
    %plot(D77(2, :), D77(1, :), 'k.--', 'MarkerSize', 23, 'LineWidth', 2)
    
end
int = interp1(ED5(:, 1),  ED5(:, 2), linspace(ED5(1,1), ED5(end,1), 100), 'spline');
if ED == 0
    plot(ED1(:, 1), ED1(:, 2),'.-', 'MarkerSize', 20, 'color', [0, 0.4470, 0.7410],  'LineWidth', 2)
    plot(ED3(:, 1), ED3(:, 2),'.-', 'MarkerSize', 20, 'color', [0, 0.4470, 0.7410],  'LineWidth', 2)
    %plot(ED5(:, 1), ED5(:, 2), 'color', [0, 0.4470, 0.7410], 'LineWidth', 2)
    plot(linspace(ED5(1,1), ED5(end,1), 100), int(:),'.-', 'MarkerSize', 20, 'color', [0, 0.4470, 0.7410],  'LineWidth', 2)
    plot(ED7(1:15, 1)-0.05, ED7(1:15, 2),'.-', 'MarkerSize', 20, 'color', [0, 0.4470, 0.7410],  'LineWidth', 2)
    %plot(ED7B(1, :), ED7B(2, :), 'color', [0, 0.4470, 0.7410], 'LineWidth', 2)
end
int2 = interp1(-DMRG5(1:end-2,3) + 0.2,  DMRG5(1:end-2,end), linspace(-DMRG5(1,3) + 0.2, -DMRG5(end-2,3) + 0.2, 50), 'spline');
if DMRG == 0
    plot(-DMRG1(1:end-3, 1), DMRG1(1:end-3, 2), 's', 'LineWidth', 1.5, 'MarkerSize', 7, 'color', [0.8500, 0.3250, 0.0980]	)
    plot(-DMRG3.alpha_list(35:2:end-16), DMRG3.delta_E_DMRG(35:2:end-16), 's', 'LineWidth', 1.5, 'MarkerSize', 7, 'color', [0.8500, 0.3250, 0.0980]	)
    %plot(-DMRG5(1:end-2,3) + 0.2, DMRG5(1:end-2,end), 's', 'LineWidth', 1.5, 'MarkerSize', 7, 'color', [0.6350, 0.0780, 0.1840])
    plot(linspace(-DMRG5(1,3) + 0.2, -DMRG5(end-2,3) + 0.2, 50), int2, 's', 'LineWidth', 1.5, 'MarkerSize', 7, 'color', [0.8500, 0.3250, 0.0980]	)
    plot(DMRG7.DMRG7Particle(1, :) + 0.1, DMRG7.DMRG7Particle(2, :), 's', 'LineWidth', 1.5, 'MarkerSize', 7, 'color', [0.8500, 0.3250, 0.0980], 'DisplayName', 'DMRG')
end

box on

set(gca, 'Yscale', 'log')

xlabel('\alpha', 'FontSize', 20)
ylabel('\Delta', 'FontSize', 20)
xlim([0 17])
ylim([0 60])
hold off


    leg{1}  = sprintf('Perpendicular Instanton factors 1 e');
    leg{2}  = sprintf('3e');
    leg{3}  = sprintf('5e');
    leg{4}  = sprintf('7e');
%     leg{5}  =sprintf('ED');
%     leg{6}  = sprintf('');
%     leg{7}  = sprintf('');
%     leg{8}  = sprintf('');
%     leg{9}  = sprintf('DMRG');
%     leg{10}  = sprintf('');
%     leg{11}  = sprintf('');
%     leg{12}  = sprintf('');


    hl=legend(leg, 'Location','NorthEast');
    legend('boxoff');
    set(hl,'Interpreter','latex', 'FontSize', 16)
    legend('boxoff');

if 1
        %set(gcf, 'PaperPositionMode', 'auto');
        fname = sprintf('Perp_fact.pdf');

        set(gcf,'paperunits','in');
        set(gcf,'papersize',[12,6]); % Desired outer dimensions
        % of figure

        hfig = gcf;
        print(hfig,'-bestfit','-dpdf',fname);


    end
