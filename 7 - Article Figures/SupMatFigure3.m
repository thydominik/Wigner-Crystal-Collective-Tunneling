clc
clear all

F = load('Freqs.mat'); F = F.Frequencies;
A = load('AlphaVals.mat'); A = A.AlphaValues;

figure(15)
clf(figure(15))
hold on
plot(A, sqrt(F), '-', 'LineWidth', 3)
axis square
box
grid
ax = gca;
ax.FontSize = 20;
xlabel('\alpha', 'FontSize', 25)
ylabel('\omega_i', 'FontSize', 25)
xline(4.443, 'LineWidth', 2)
hold off