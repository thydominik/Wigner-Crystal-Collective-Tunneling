clc
clear all

n1_d = load('n1.csv');
n3_d = load('n3.csv');
n5_d = load('n5.csv');
n7_d = load('n7.csv');

figure(1)
clf(figure(1))
hold on
%title('experimental data')
xlabel('V_{k2} [mV]', 'FontSize', 12)
ylabel('\Delta [a.u.]', 'FontSize', 12)
plot(n1_d(:, 1), n1_d(:, 2), '-', 'LineWidth', 2, 'DisplayName', '1 particle')
plot(n3_d(:, 1), n3_d(:, 2), '-', 'LineWidth', 2, 'DisplayName', '3 particle')
plot(n5_d(:, 1), n5_d(:, 2), '-', 'LineWidth', 2, 'DisplayName', '5 particle')
plot(n7_d(:, 1), n7_d(:, 2), '-', 'LineWidth', 2, 'DisplayName', '7 particle')
legend
hold off

%scaling:
ED1 = load("splitting energy from Schrödinger.mat");
ED1 = ED1.dE;
%%
param1_a = 18;
param1_d = 30.5;
figure(2)
clf(figure(2))
hold on
%title('1 part. exp & ED')
xlabel('V_{k2} or \alpha_{scaled}', 'FontSize', 15)
ylabel('\Delta_{scaled} or t[a.u.]', 'FontSize', 15)
plot(n1_d(:, 1), n1_d(:, 2), 'o-', 'DisplayName', 'Experimental data')
plot((ED1(:, 1)) * param1_a + 68, ED1(:, 2) * param1_d, '.-', 'DisplayName', 'Schrödinger ED')
text(150,20,'\alpha to V_{k2} => \alpha * 18 + 68 = V_{k2}')
text(150,15,'\Delta to exp t[a.u.] => \Delta * 30.5 = t')
legend
hold off
%%
figure(3)
clf(figure(3))
hold on
title('Scaling for 1 particle')
xlabel('\alpha', 'FontSize', 15)
ylabel('V_{k2}', 'FontSize', 15)
plot(ED1(:, 1), ED1(:, 1) * param1_a, 'o-')
hold off
%%
%n3_d(:, 1) = n3_d(:, 1) - min(n3_d(:, 1));
ED3     = load('E_Schrodinger_3e_eta_20.00_beta_0.01_N_100.dat');
dDE3    = abs(ED3(:, 3) - ED3(:, 2));
alpha   = -ED3(:, 1);
param3_a    = 95;
param3_d    = 52;
crit_alpha  = 4.23;
figure(4)
clf(figure(4))
hold on
title('3 partical experimental data & ED')
xlabel('V_{k2} or \alpha_{scaled}', 'FontSize', 15)
ylabel('\Delta_{scaled} or t[a.u.]', 'FontSize', 15)
plot(n3_d(:, 1), n3_d(:, 2), 'o-', 'DisplayName', 'Experimental data')
%From eq. pos the alpha critical is 4.23 ish
plot((alpha - crit_alpha) * param3_a , dDE3 * param3_d, '.-', 'DisplayName', 'Schrödinger ED')
text(350,30,'\alpha to V_{k2} => (\alpha - \alpha_c) * 95 = V_{k2}')
text(350,20,'\Delta to exp t[a.u.] => \Delta * 52 = t')
ylim([0 50])
legend
hold off
%%
figure(5)
clf(figure(5))
hold on
title('Scaling for 3 particles')
xlabel('\alpha - min(\alpha)', 'FontSize', 15)
ylabel('V_{k2}', 'FontSize', 15)
plot(alpha - crit_alpha, (alpha - crit_alpha) * param3_a, 'o-')
hold off
%%
figure(6)
clf(figure(6))
hold on
%title('a = 1.16; b = 0.9')
plot(ED1(:, 1), ED1(:, 2), '-', 'LineWidth', 2, 'DisplayName', '1 part. - ED')
plot((alpha), dDE3 , '-', 'LineWidth', 2, 'DisplayName', '3 part. - ED')
%set(gca, 'Yscale', 'log')
plot((n1_d(:, 1) -68) * (1/param1_a), n1_d(:, 2) * (1/param1_d), 'o-', 'LineWidth', 1.2, 'DisplayName', '1 part. - Exp')
plot((n3_d(:, 1 ) + 410) * (1/param3_a), n3_d(:, 2) * (1/param3_d), 'o-', 'LineWidth', 1.2, 'DisplayName', '3 part. - Exp')
xlabel('\alpha ', 'FontSize', 15)
ylabel('\Delta', 'FontSize', 15)

xlim([0, 10])
legend
hold off

%%
figure(7)
clf(figure(7))
hold on
title('experimental data rescaled')
xlabel('V_{k2} [mV]')
ylabel('t [a.u.]')
V = [110 330 610 850];
V = V;
plot(n1_d(:, 1) - V(1)     , n1_d(:, 2) * 1, 'o')
plot((n3_d(end/2:end, 1) - V(2)) * 0.8  , n3_d(end/2:end, 2) * 6, 'x')
plot((n5_d(:, 1) - V(3))* 0.55  , n5_d(:, 2) * 9, '.')
plot((n7_d(:, 1) - V(4))* 0.46 , n7_d(:, 2) * 6.5, '*')
%xlim([40 1200])
%xlim([40 160])
%ylim([0 50])
hold off

%%
Inst1 = load('splitting energy from Instanton calc.mat');
Inst1 = Inst1.splitt;

figure(8)
clf(figure(8))
hold on
%plot(0.995 * abs(Inst1(:, 1)) * param1_a + 68,  Inst1(:, 2) * param1_d , '.-', 'DisplayName', 'Instanton same fit')
plot(abs(Inst1(:, 1)) * 20 + 70,  Inst1(:, 2) * 15 , 'b.-', 'MarkerSize', 15, 'DisplayName', 'Instanton new fit')
plot(abs(ED1(:, 1)) * param1_a + 68, abs(ED1(:, 2)) * param1_d, 'r-','LineWidth', 2, 'DisplayName', 'ED')
plot(n1_d(:, 1), n1_d(:, 2), 'ko-', 'MarkerFaceColor', 'k', 'LineWidth', 1, 'DisplayName', 'Experimental data')
set(gca, 'Yscale', 'log')
legend
grid on
ylim([0 35])
xlim([50 300])
xlabel('V_{k2}', 'FontSize', 15)
ylabel('\Delta', 'FontSize', 15)
hold off

%%
Inst1 = load('splitting energy from Instanton calc.mat');
Inst1 = Inst1.splitt;

figure(9)
clf(figure(9))
hold on
plot(0.995 * abs(Inst1(:, 1)) * param1_a + 68,  Inst1(:, 2) * param1_d , 'b.-', 'MarkerSize', 15, 'DisplayName', 'Instanton same fit')
%plot(abs(Inst1(:, 1)) * 20 + 70,  Inst1(:, 2) * 15 , '.-', 'DisplayName', 'Instanton new fit')
plot(abs(ED1(:, 1)) * param1_a + 68, abs(ED1(:, 2)) * param1_d, 'r-','LineWidth', 2, 'DisplayName', 'ED')
plot(n1_d(:, 1), n1_d(:, 2), 'ko-', 'MarkerFaceColor', 'k', 'LineWidth', 1, 'DisplayName', 'Experimental data')
set(gca, 'Yscale', 'log')
legend
grid on
ylim([0 35])
xlim([50 300])
xlabel('V_{k2}', 'FontSize', 15)
ylabel('\Delta', 'FontSize', 15)
hold off

%%
Inst3 = load('splitting energy from Instanton calc 3 part.mat');
Inst3 = Inst3.data;

figure(10)
clf(figure(10))
hold on
title('3 partical experimental data & ED')
xlabel('V_{k2} or \alpha_{scaled}', 'FontSize', 15)
ylabel('\Delta_{scaled} or t[a.u.]', 'FontSize', 15)
plot(n3_d(:, 1), n3_d(:, 2), 'ko-', 'MarkerFaceColor', 'k', 'LineWidth', 1, 'DisplayName', 'Experimental data')
plot( (abs(Inst3(:, 1)) - crit_alpha) * param3_a,  param3_d * Inst3(:, 2) .* Inst3(:, 8),'b.-', 'MarkerSize', 15, 'DisplayName', 'Instanton same fit')
%plot( (abs(Inst3(:, 1)) - crit_alpha) * 100, 25 * Inst3(:, 2) .* Inst3(:, 8),'o-', 'DisplayName', 'Instanton new fit')
%From eq. pos the alpha critical is 4.23 ish
plot((alpha - crit_alpha) * param3_a , dDE3 * param3_d, 'r-', 'LineWidth', 2, 'DisplayName', 'ED')
%text(350,30,'\alpha to V_{k2} => (\alpha - \alpha_c) * 95 = V_{k2}')
%text(350,20,'\Delta to exp t[a.u.] => \Delta * 52 = t')
ylim([0 50])
legend
%set(gca, 'Yscale', 'log')
ylim([0 60])
xlim([0 700])
xlabel('V_{k2}', 'FontSize', 15)
ylabel('\Delta', 'FontSize', 15)
hold off
%%
Inst3 = load('splitting energy from Instanton calc 3 part.mat');
Inst3 = Inst3.data;

figure(11)
clf(figure(11))
hold on
title('3 partical experimental data & ED')
xlabel('V_{k2} or \alpha_{scaled}', 'FontSize', 15)
ylabel('\Delta_{scaled} or t[a.u.]', 'FontSize', 15)
plot(n3_d(:, 1), n3_d(:, 2), 'ko-', 'MarkerFaceColor', 'k', 'LineWidth', 1, 'DisplayName', 'Experimental data')
%plot( (abs(Inst3(:, 1)) - crit_alpha) * param3_a, 0.6 * param3_d * Inst3(:, 2) .* Inst3(:, 8),'.-', 'DisplayName', 'Instanton same fit')
plot( (abs(Inst3(:, 1)) - crit_alpha) * 100, 25 * Inst3(:, 2) .* Inst3(:, 8),'b.-', 'MarkerSize', 15, 'DisplayName', 'Instanton new fit')
%From eq. pos the alpha critical is 4.23 ish
plot((alpha - crit_alpha) * param3_a , dDE3 * param3_d, 'r-', 'LineWidth', 2, 'DisplayName', 'ED')
%text(350,30,'\alpha to V_{k2} => (\alpha - \alpha_c) * 95 = V_{k2}')
%text(350,20,'\Delta to exp t[a.u.] => \Delta * 52 = t')
ylim([0 50])
legend
set(gca, 'Yscale', 'log')
ylim([0 60])
xlim([0 700])
xlabel('V_{k2}', 'FontSize', 15)
ylabel('\Delta', 'FontSize', 15)
hold off
