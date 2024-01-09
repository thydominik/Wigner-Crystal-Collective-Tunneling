%% eta = 4
close all
Ne      = 3;
Nvalues = 20;
Norb    = 8;
Norb_s  = 8;
eta     = 2;
Sz1     = 1.5;
Sz2     = 0.5;
levels  = 6;

Et = load(sprintf('Energy_spectrum_Sz_%.1f_Ne_%d_Norb_%d_eta_%.4f.dat', Sz1, Ne, Norb, eta));
Es = load(sprintf('Energy_spectrum_Sz_%.1f_Ne_%d_Norb_%d_eta_%.4f.dat', Sz2, Ne, Norb_s, eta));

E_sch = load('E_Schrodinger_3e_spinless_eta_20.00_N_200_beta_0.30_eigs.dat');

h = figure('color', 'white', 'units', 'inches', 'position', [1 1 8 6]);

hold on
colors      = ['r', 'k', 'b', 'm', 'y', 'g'];
colors2     = ['r', 'r', 'r', 'k', 'k', 'k'];

fontsize(h, 12, "points")

Emin        = min(Es(:, 2), Et(:, 2));

symbols     = ['o', 's', 'd', '^', 'v', 'x'];
symbolsT    = ['x', 's', 'd', '^', 'v', 'x'];

plot(-Es(:, 1), Es(:, 6 + 1) - Emin(:), symbols(6), 'LineWidth', 1, 'MarkerEdgeColor', colors2(6), 'MarkerSize', 6);
plot(-Es(:, 1), Es(:, 5 + 1) - Emin(:), symbols(5), 'LineWidth', 1, 'MarkerEdgeColor', colors2(5), 'MarkerSize', 6);
plot(-Es(:, 1), Es(:, 4 + 1) - Emin(:), symbols(4), 'LineWidth', 1, 'MarkerEdgeColor', colors2(4), 'MarkerSize', 6);
plot(-Es(:, 1), Es(:, 3 + 1) - Emin(:), symbols(3), 'LineWidth', 1, 'MarkerEdgeColor', colors2(3), 'MarkerSize', 6);
plot(-Es(:, 1), Es(:, 2 + 1) - Emin(:), symbols(2), 'LineWidth', 1, 'MarkerEdgeColor', colors2(2), 'MarkerSize', 6);
plot(-Es(:, 1), Es(:, 1 + 1) - Emin(:), symbols(1), 'LineWidth', 1, 'MarkerEdgeColor', colors2(1), 'MarkerSize', 6);

xlabel(sprintf('$\\alpha$'),    'Interpreter', 'latex', 'FontSize', 30);
ylabel(sprintf('$E_n$'),        'Interpreter', 'latex', 'FontSize', 30);

xlim([2, 9])
ylim([0 2.5])

box on
axis square

legend 'DMRG - S_{TOT} = 3/2' 'DMRG - S_{TOT} = 3/2' 'DMRG - S_{TOT} = 3/2' 'DMRG - S_{TOT} = 1/2' 'DMRG - S_{TOT} = 1/2' 'DMRG - S_{TOT} = 1/2'
text(2.5, 0.5, '$\eta = 4$', 'Interpreter', 'latex', FontSize=18)
if 1
    set(gcf, 'PaperPositionMode', 'auto');
    fname = sprintf('Spinless_Energy_Ne_%d_Norb_%d_%d_eta_%.2f.pdf', Ne, Norb * Ne, eta);

    set(gcf, 'paperunits', 'in');
    set(gcf, 'papersize', [8.2,6.4]); % Desired outer dimensions of figure

    hfig = gcf;

    print(hfig, '-dpdf', fname);
end

%% eta = 20

close all
Ne      = 3;
Nvalues = 20;
Norb    = 8;
Norb_s  = 8;
eta     = 10;
Sz1     = 1.5;
Sz2     = 0.5;
levels  = 6;

Et = load(sprintf('Energy_spectrum_Sz_%.1f_Ne_%d_Norb_%d_eta_%.4f.dat', Sz1, Ne, Norb, eta));
Es = load(sprintf('Energy_spectrum_Sz_%.1f_Ne_%d_Norb_%d_eta_%.4f.dat', Sz2, Ne, Norb_s, eta));

E_sch = load('E_Schrodinger_3e_spinless_eta_20.00_N_200_beta_0.30_eigs.dat');

h = figure('color', 'white', 'units', 'inches', 'position', [1 1 8 6]);

hold on
colors      = ['r', 'k', 'b', 'm', 'y', 'g'];
colors2     = ['r', 'r', 'r', 'k', 'k', 'k'];

fontsize(h, 12, "points")

Emin        = min(Es(:, 2), Et(:, 2));

symbols     = ['o', 's', 'd', '^', 'v', 'x'];
symbolsT    = ['x', 's', 'd', '^', 'v', 'x'];

plot(-Et(:, 1), Et(:, 2 + 1) - Emin(:), symbolsT(1), 'LineWidth', 1, 'MarkerEdgeColor', colors(2), 'MarkerSize', 6);
plot(-Es(:, 1), Es(:, 5 + 1) - Emin(:), symbols(5), 'LineWidth', 1, 'MarkerEdgeColor', colors2(5), 'MarkerSize', 6);
plot(-Es(:, 1), Es(:, 1 + 1) - Emin(:), symbols(1), 'LineWidth', 1, 'MarkerEdgeColor', colors2(1), 'MarkerSize', 6);
plot(-Es(:, 1), Es(:, 2 + 1) - Emin(:), symbols(2), 'LineWidth', 1, 'MarkerEdgeColor', colors2(2), 'MarkerSize', 6);
plot(-Et(:, 1), Et(:, 1 + 1) - Emin(:), symbols(4), 'LineWidth', 1, 'MarkerEdgeColor', colors(1), 'MarkerSize', 6);


plot(-E_sch(1:1:end, 1), E_sch(1:1:end, 8) - E_sch(1:1:end, 2), 'b-', 'LineWidth', 1.5 , 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'MarkerSize', 4);

%set(narrow_colorbar,'FontSize',fsize); %adds a scaling bar

xlabel(sprintf('$\\alpha$'),    'Interpreter', 'latex', 'FontSize', 30);
ylabel(sprintf('$E_n$'),        'Interpreter', 'latex', 'FontSize', 30);
xlim([2, 9]);

ylim([0 2.5])
box on
axis square
legend 'DMRG - S_{TOT} = 3/2' 'DMRG - S_{TOT} = 1/2' 'DMRG - S_{TOT} = 1/2' 'DMRG - S_{TOT} = 1/2' 'DMRG - S_{TOT} = 3/2' 'ED - Spinless'
text(2.5, 0.5, '$\eta = 20$', 'Interpreter', 'latex', FontSize=18)
if 1
    set(gcf, 'PaperPositionMode', 'auto');
    fname = sprintf('Spinless_Energy_Ne_%d_Norb_%d_%d_eta_%.2f.pdf', Ne, Norb * Ne, eta);

    set(gcf, 'paperunits', 'in');
    set(gcf, 'papersize', [8.2,6.4]); % Desired outer dimensions of figure

    hfig = gcf;

    print(hfig, '-dpdf', fname);
end