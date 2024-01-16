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
E_scale = 5.57;
h = figure('color', 'white', 'units', 'inches', 'position', [1 1 7 5]);

hold on
colors      = ['r', 'k', 'b', 'm', 'y', 'g'];
colors2     = ['r', 'r', 'k', 'r', 'r', 'k'];

fontsize(h, 20, "points")

Emin        = min(Es(:, 2), Et(:, 2));

symbols     = ['o', 's', 'd', '^', 'v', 'x'];
symbolsT    = ['x', 's', 'd', '^', 'v', 'x'];

plot(-Es(:, 1), E_scale*(Es(:, 6 + 1) - Emin(:)), symbols(6), 'LineWidth', 1, 'MarkerEdgeColor', colors2(6), 'MarkerSize', 6);
plot(-Es(:, 1), E_scale*(Es(:, 5 + 1) - Emin(:)), symbols(5), 'LineWidth', 1, 'MarkerEdgeColor', colors2(5), 'MarkerSize', 6);
plot(-Es(:, 1), E_scale*(Es(:, 4 + 1) - Emin(:)), symbols(4), 'LineWidth', 1, 'MarkerEdgeColor', colors2(4), 'MarkerSize', 6);
plot(-Es(:, 1), E_scale*(Es(:, 3 + 1) - Emin(:)), symbols(3), 'LineWidth', 1, 'MarkerEdgeColor', colors2(3), 'MarkerSize', 6);
plot(-Es(:, 1), E_scale*(Es(:, 2 + 1) - Emin(:)), symbols(2), 'LineWidth', 1, 'MarkerEdgeColor', colors2(2), 'MarkerSize', 6);
plot(-Es(:, 1), E_scale*(Es(:, 1 + 1) - Emin(:)), symbols(1), 'LineWidth', 1, 'MarkerEdgeColor', colors2(1), 'MarkerSize', 6);

xlabel(sprintf('$\\alpha$'),    'Interpreter', 'latex', 'FontSize', 30);
ylabel(sprintf('$E_n (K)$'),        'Interpreter', 'latex', 'FontSize', 30);

xlim([2, 9])
ylim([1e-3 2.5*E_scale])
%set(gca, 'YScale', 'log');
box on
%axis square

leg=legend ('DMRG - S_{TOT} = 3/2', 'DMRG - S_{TOT} = 1/2' ,'DMRG - S_{TOT} = 1/2' ,'DMRG - S_{TOT} = 3/2' ,'DMRG - S_{TOT} = 1/2' ,'DMRG - S_{TOT} = 1/2', 'FontSize', 16)
leg.Position =  [0.3393 0.5125 0.3611 0.4681];
text(2.5, 0.5*E_scale, '$\eta = 4$', 'Interpreter', 'latex', FontSize=18)
if 1
    set(gcf, 'PaperPositionMode', 'auto');
    fname = sprintf('Spinhalf_Energy_Ne_%d_Norb_%d_eta_%.2f.pdf', Ne, Norb * Ne, 2*eta);

    set(gcf, 'paperunits', 'in');
    set(gcf, 'papersize', [7.5,5.5]); % Desired outer dimensions of figure

    hfig = gcf;

    print(hfig, '-dpdf', fname);
end
%% eta = 20

%close all
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

h = figure('color', 'white', 'units', 'inches', 'position', [1 1 7 5]);

hold on
colors      = ['r', 'k', 'b', 'm', 'y', 'g'];
colors2     = ['r', 'r', 'r', 'k', 'k', 'k'];

fontsize(h, 20, "points")

Emin        = min(Es(:, 2), Et(:, 2));

symbols     = ['o', 's', 'd', '^', 'v', 'x'];
symbolsT    = ['x', 's', 'd', '^', 'v', 'x'];

plot(-Et(:, 1), E_scale*(Et(:, 2 + 1) - Emin(:)), symbolsT(1), 'LineWidth', 1, 'MarkerEdgeColor', colors(2), 'MarkerSize', 6);

plot(-Es(:, 1), E_scale*(Es(:, 5 + 1) - Emin(:)), symbols(5), 'LineWidth', 1, 'MarkerEdgeColor', colors2(5), 'MarkerSize', 6);
plot(-Es(:, 1), E_scale*(Es(:, 6) - Emin(:)), symbols(4), 'LineWidth', 1, 'MarkerEdgeColor', colors2(5), 'MarkerSize', 6);
plot(-Es(:, 1), E_scale*(Es(:, 1 + 1) - Emin(:)), symbols(1), 'LineWidth', 1, 'MarkerEdgeColor', colors2(1), 'MarkerSize', 6);
plot(-Es(:, 1), E_scale*(Es(:, 2 + 1) - Emin(:)), symbols(2), 'LineWidth', 1, 'MarkerEdgeColor', colors2(2), 'MarkerSize', 6);
plot(-Et(:, 1), E_scale*(Et(:, 1 + 1) - Emin(:)), symbols(4), 'LineWidth', 1, 'MarkerEdgeColor', colors(1), 'MarkerSize', 6);



plot(-E_sch(1:1:end, 1), E_scale*(E_sch(1:1:end, 8) - E_sch(1:1:end, 2)), 'b-', 'LineWidth', 1.5 , 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'MarkerSize', 4);

%set(narrow_colorbar,'FontSize',fsize); %adds a scaling bar

xlabel(sprintf('$\\alpha$'),    'Interpreter', 'latex', 'FontSize', 30);
ylabel(sprintf('$E_n(K)$'),        'Interpreter', 'latex', 'FontSize', 30);
xlim([2, 9]);
%set(gca, 'YScale', 'log');
ylim([1e-6 4.0*E_scale]);

box on
%axis square
legend('DMRG - S_{TOT} = 3/2', 'DMRG - S_{TOT} = 1/2', 'DMRG - S_{TOT} = 1/2', 'DMRG - S_{TOT} = 1/2' ,'DMRG - S_{TOT} = 1/2' ,'DMRG - S_{TOT} = 3/2', 'ED - Spinless','Fontsize', 16);
text(2.5, 0.5*E_scale, '$\eta = 20$', 'Interpreter', 'latex', FontSize=18)
if 1
    set(gcf, 'PaperPositionMode', 'auto');
    fname = sprintf('Spinhalf_Energy_Ne_%d_Norb_%d_eta_%.2f.pdf', Ne, Norb * Ne, 2*eta);

    set(gcf, 'paperunits', 'in');
    set(gcf, 'papersize', [7.5,5.5]); % Desired outer dimensions of figure

    hfig = gcf;

    print(hfig, '-dpdf', fname);
end