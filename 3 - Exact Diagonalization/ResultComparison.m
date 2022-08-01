clc
clear all

disp('1 particles')


a1 = load('EDSplitting_1_particle_Nx_50.mat');
a2 = load('EDSplitting_1_particle_Nx_100.mat');
a3 = load('EDSplitting_1_particle_Nx_200.mat');
a4 = load('EDSplitting_1_particle_Nx_400.mat');
a5 = load('EDSplitting_1_particle_Nx_800.mat');
dmrg = load('OneParticleDMRG.txt');

figure(1)
clf(figure(1))
hold on
title('')
xlabel('\alpha')
ylabel('\Delta')
legend
%plot(a1.data(:, 1), a1.data(:, 2), '.-', 'DisplayName', 'Nx = 10')
%plot(a2.data(:, 1), a2.data(:, 2), '.-', 'DisplayName', 'Nx = 20')
plot(a1.data(:, 1), a1.data(:, 2), '.-', 'DisplayName', 'Nx = 50')
%plot(a2.data(:, 1), a2.data(:, 2), '.-', 'DisplayName', 'Nx = 100')
%plot(a3.data(:, 1), a3.data(:, 2), '.-', 'DisplayName', 'Nx = 200')
%plot(a4.data(:, 1), a4.data(:, 2), '.-', 'DisplayName', 'Nx = 400')
plot(a5.data(:, 1), a5.data(:, 2), '.-', 'DisplayName', 'Nx = 800')
plot(-dmrg(:, 1), dmrg(:, 2), 'o', 'DisplayName', 'DMRG')
set(gca, 'YScale', 'log')
hold off

%%

disp('2 particles')

a1 = load('EDSplitting_2_particles_Nx_10.mat');
a2 = load('EDSplitting_2_particles_Nx_20.mat');
a3 = load('EDSplitting_2_particles_Nx_50.mat');
a4 = load('EDSplitting_2_particles_Nx_100.mat');
a5 = load('EDSplitting_2_particles_Nx_200.mat');
a6 = load('EDSplitting_2_particles_Nx_300.mat');
a1 = load('EDSplitting_2_particles_restricted_Nx1_30_Nx2_30_beta_1e-05.mat');
a2 = load('EDSplitting_2_particles_restricted_Nx1_50_Nx2_50_beta_1e-05.mat');
a3 = load("EDSplitting_2_particles_restricted_Nx1_70_Nx2_70_beta_1e-05.mat");
a4 = load('EDSplitting_2_particles_restricted_Nx1_100_Nx2_100_beta_1e-05.mat');
a5 = load("EDSplitting_2_particles_restricted_Nx1_200_Nx2_200_beta_1e-05.mat");

figure(1)
clf(figure(1))
hold on
title('')
xlabel('\alpha')
ylabel('\Delta')
legend
%plot(a1.data.alpha, a1.data.EnergySplitting, '.-', 'DisplayName', 'Nx1 = Nx2 = 30')
%%plot(a2.data.alpha, a2.data.EnergySplitting, '.-', 'DisplayName', 'Nx1 = Nx2 = 50')
%%plot(a3.data.alpha, a3.data.EnergySplitting, '.-', 'DisplayName', 'Nx1 = Nx2 = 70')
%plot(a4.data.alpha, a4.data.EnergySplitting, '.-', 'DisplayName', 'Nx1 = Nx2 = 100')
plot(a5.data.alpha, a5.data.EnergySplitting, '.-', 'DisplayName', 'Nx1 = Nx2 = 200')
plot(a6.data(:, 1), a6.data(:, 2) - a5.data.EnergySplitting, 'o-', 'DisplayName', 'Nx = 300')
set(gca, 'YScale', 'log')
hold off

%%
disp('2 particles')

a1 = load('EDSplitting_3_particles_Nx_10.mat');
a2 = load('EDSplitting_3_particles_Nx_20.mat');
a3 = load('EDSplitting_3_particles_Nx_30.mat');
a4 = load('EDSplitting_3_particles_Nx_40.mat');


figure(1)
clf(figure(1))
hold on
plot(a1.data(:, 1), a1.data(:, 2), '.-', 'DisplayName', 'Nx = 10')
plot(a2.data(:, 1), a2.data(:, 2), '.-', 'DisplayName', 'Nx = 20')
plot(a3.data(:, 1), a3.data(:, 2), '.-', 'DisplayName', 'Nx = 30')
plot(a4.data(:, 1), a4.data(:, 2), '.-', 'DisplayName', 'Nx = 40')
title('')
xlabel('\alpha')
ylabel('\Delta')
legend
set(gca, 'YScale', 'log')
hold off