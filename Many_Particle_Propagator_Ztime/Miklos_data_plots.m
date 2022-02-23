%Plotting

M_data1 = load('E_Schrodinger_3e_eta_20.00_beta_0.01_N_100.dat');
M_data2 = load('E_Schrodinger_3e_eta_20.00_beta_0.10_N_100.dat');
M_data3 = load('E_Schrodinger_3e_eta_20.00_beta_0.30_N_60.dat');
M_data4 = load('E_Schrodinger_3e_eta_20.00_beta_0.60_N_100.dat');

figure(5)
clf(figure(5))
hold on
plot(-M_data1(:, 1), M_data1(:, 3) - M_data1(:, 2), '.-', 'DisplayName', 'Beta = 0.01')
plot(-M_data2(:, 1), M_data2(:, 3) - M_data2(:, 2), '.-', 'DisplayName', 'Beta = 0.10')
plot(-M_data3(:, 1), M_data3(:, 3) - M_data3(:, 2), '.-', 'DisplayName', 'Beta = 0.30')
plot(-M_data4(:, 1), M_data4(:, 3) - M_data4(:, 2), '.-', 'DisplayName', 'Beta = 0.60')
xlim([5 10])
ylim([10^-5.2 2])
set(gca, 'Yscale', 'log')
xlabel('- \alpha', 'FontSize', 24)
ylabel('\Delta_{ED}', 'FontSize', 24)
legend
hold off