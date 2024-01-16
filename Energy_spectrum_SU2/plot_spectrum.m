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

h = figure('color', 'white', 'units', 'inches', 'position', [1 1 10 10]);

hold on
colors      = ['r', 'k', 'b', 'm', 'y', 'g'];
colors2     = ['r', 'r', 'r', 'k', 'k', 'k'];

fontsize(h, 12, "points")

Emin        = min(Es(:, 2), Et(:, 2));

symbols     = ['o', 's', 'd', '^', 'v', 'x'];
symbolsT    = ['x', 's', 'd', '^', 'v', 'x'];

for j = 1:5%[1, 2,  5]

    plot(-Es(:, 1), Es(:, j + 1) - Emin(:), symbols(j), 'LineWidth', 1, 'MarkerEdgeColor', colors2(j), 'MarkerSize', 6);
    %plot(-Es(:, 1), Es(:, j + 1) - Emin(:), symbols(j), 'LineWidth',  1, 'MarkerEdgeColor', colors2(j),  'MarkerSize', 6);
end

for j = 1%:levels-4
    %plot(-Et(:, 1), Et(:, j + 1) - Emin(:), symbolsT(j), 'LineWidth', 1, 'MarkerEdgeColor', colors(j), 'MarkerSize', 6);
end

%plot(-E_sch(1:2:end, 1), E_sch(1:2:end, 8) - E_sch(1:2:end, 2), 'b-', 'LineWidth', 1, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'MarkerSize', 4);

%set(narrow_colorbar,'FontSize',fsize); %adds a scaling bar

xlabel(sprintf('$\\alpha$'),    'Interpreter', 'latex', 'FontSize', 30);
ylabel(sprintf('$E_n$'),        'Interpreter', 'latex', 'FontSize', 30);
xlim([2, 9]);

ylim([0 2.5])
box on
axis square
legend 'DMRG - S_{TOT} = 1/2' 'DMRG - S_{TOT} = 1/2' 'DMRG - S_{TOT} = 3/2' 'DMRG - S_{TOT} = 3/2' 'DMRG - S_{TOT} = 1/2' 'DMRG - S_{TOT} = 3/2' 'ED - Spinless'
text(2.5, 0.5, '$\eta = 20$', 'Interpreter', 'latex', FontSize=18)
if 1
    set(gcf, 'PaperPositionMode', 'auto');
    fname = sprintf('Spinless_Energy_Ne_%d_Norb_%d_%d_eta_%.2f.pdf', Ne, Norb * Ne, eta);

    set(gcf, 'paperunits', 'in');
    set(gcf, 'papersize', [8.2,6.4]); % Desired outer dimensions of figure

    hfig = gcf;

    print(hfig, '-dpdf', fname);
end
