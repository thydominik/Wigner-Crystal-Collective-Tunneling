clc
clear all

addpath('D:\BME PhD\Wigner Crystal Collective Tunneling\Data\Equilibrium Positions')

d1 = load("Eq_Pos_eta_20_particles_1.mat");
d1 = d1.eqpos;

d3 = load("Eq_Pos_eta_20_particles_3.mat");
d3 = d3.eqpos;

d5 = load("Eq_Pos_eta_20_particles_5.mat");
d5 = d5.eqpos;

d7 = load("Eq_Pos_eta_20_particles_7.mat");
d7 = d7.eqpos;

figure(1)
clf(figure(1))
hold on
plot(d1(2, :), d1(1, :), '.')
xline(0)
hold off

figure(2)
clf(figure(2))
hold on
plot(d3(end, :), d3(1, :), '.')
plot(d3(end, :), d3(2, :), '.')
plot(d3(end, :), d3(3, :), '.')
xline(-4.45)
hold off

figure(3)
clf(figure(3))
hold on
plot(d5(end, :), d5(1, :), '.')
plot(d5(end, :), d5(2, :), '.')
plot(d5(end, :), d5(3, :), '.')
plot(d5(end, :), d5(4, :), '.')
plot(d5(end, :), d5(5, :), '.')
xline(-7.81)
hold off

figure(4)
clf(figure(4))
hold on
plot(d7(end, :), d7(1, :), '.')
plot(d7(end, :), d7(2, :), '.')
plot(d7(end, :), d7(3, :), '.')
plot(d7(end, :), d7(4, :), '.')
plot(d7(end, :), d7(5, :), '.')
plot(d7(end, :), d7(6, :), '.')
plot(d7(end, :), d7(7, :), '.')
xline(-10.61)
hold off

figure(5)
clf(figure(5))
hold on
ac = [0 -4.45 -7.81 -10.61 -13.0625 -15.2];
pn = [1 3 5 7 9 11];
x = linspace(1, 11, 100);
func = 5.151 - 5.131 * x.^(0.5756);
plot(pn, -ac, 'o', 'DisplayName', 'MC classical critical values')
plot(x, -func, '.', 'DisplayName', 'fit of form a+b*x^c')
xlabel('Particle number')
ylabel('|\alpha_c|')
grid on
legend
xlim([1 11])
ylim([0 16])
hold off
