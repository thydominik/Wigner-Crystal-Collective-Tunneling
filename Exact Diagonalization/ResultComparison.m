clc
clear all

disp('1 particles')

a1 = load('EDSplitting_1_particle_Nx_10.mat');
a2 = load('EDSplitting_1_particle_Nx_20.mat');
a3 = load('EDSplitting_1_particle_Nx_50.mat');
a4 = load('EDSplitting_1_particle_Nx_100.mat');
a5 = load('EDSplitting_1_particle_Nx_200.mat');
a6 = load('EDSplitting_1_particle_Nx_500.mat');
a7 = load('EDSplitting_1_particle_Nx_1000.mat');
figure(1)
clf(figure(1))
hold on
title('')
xlabel('\alpha')
ylabel('\Delta')
legend
%plot(a1.data(:, 1), a1.data(:, 2), '.-', 'DisplayName', 'Nx = 10')
%plot(a2.data(:, 1), a2.data(:, 2), '.-', 'DisplayName', 'Nx = 20')
plot(a3.data(:, 1), a3.data(:, 2), '.-', 'DisplayName', 'Nx = 50')
plot(a4.data(:, 1), a4.data(:, 2), '.-', 'DisplayName', 'Nx = 100')
plot(a5.data(:, 1), a5.data(:, 2), '.-', 'DisplayName', 'Nx = 200')
plot(a6.data(:, 1), a6.data(:, 2), '.-', 'DisplayName', 'Nx = 500')
plot(a7.data(:, 1), a7.data(:, 2), '.-', 'DisplayName', 'Nx = 1000')
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


figure(1)
clf(figure(1))
hold on
title('')
xlabel('\alpha')
ylabel('\Delta')
legend
%plot(a1.data(:, 1), a1.data(:, 2), '.-', 'DisplayName', 'Nx = 10')
%plot(a2.data(:, 1), a2.data(:, 2), '.-', 'DisplayName', 'Nx = 20')
plot(a3.data(:, 1), a3.data(:, 2), '.-', 'DisplayName', 'Nx = 50')
plot(a4.data(:, 1), a4.data(:, 2), '.-', 'DisplayName', 'Nx = 100')
plot(a5.data(:, 1), a5.data(:, 2), '.-', 'DisplayName', 'Nx = 200')
plot(a6.data(:, 1), a6.data(:, 2), '.-', 'DisplayName', 'Nx = 300')
set(gca, 'YScale', 'log')
hold off
