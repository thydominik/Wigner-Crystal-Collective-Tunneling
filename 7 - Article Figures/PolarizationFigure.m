%% Figure 4: Polarization 3D plot theoretical
% 
clear all
close all
clc

FontSize = 20;
Position = [1 1 6 5];

Figure = figure('color', 'white', 'units', 'inches', 'Position', Position);

hold on
% Fontsize:
    set(gca, 'FontSize', FontSize)

% Axis limits:
    %xlim([5 9])
    %xticks([-7.5 -5 -2.5 0 2.5 5 7.5])
    %ylim([-0.1 0.1])
    yticks([4 5 6 7 8 9 10 11 12 13 14 15])
    %zlim([])
    %axis square
box on
% Axis Labels:
title('$\tilde{\mathcal{P}}$', 'Interpreter', 'latex', 'FontSize', FontSize)
    xlabel('$\epsilon$', 'Interpreter', 'latex', 'FontSize', FontSize)
    ylabel('$\alpha$', 'FontWeight', 'Bold', 'Interpreter', 'latex', 'FontWeight', 'Bold', 'FontSize', FontSize + 5)
    zlabel('$P(\alpha, \epsilon)$', 'Interpreter', 'latex', 'FontSize', FontSize+5)
% Data: 
    Polarization = load('Polarization3Particle.mat');
    Polarization = Polarization.Polarization;
% Plotting:
hold on
box
     surf( Polarization.Kappa, Polarization.Alpha, Polarization.Polarization , 'EdgeColor', 'interp', 'FaceColor', 'interp')
colorbar
c                   = colorbar('eastoutside');
    c                   = colorbar('eastoutside'); %, 'Direction','reverse');
    %c.Label.String      = '$\tilde{\mathcal{P}}$';
    %c.Label.Interpreter = 'latex';
%text( -0.125, 4, 0, '(b)', 'HorizontalAlignment', 'center', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize);
%     for alphaInd = 1:26
%         plot3(Polarization.Alpha(alphaInd) * ones(201, 1), Polarization.Kappa, Polarization.Polarization(alphaInd, :, :), 'Color', [alphaInd/26 0.5 (1 - alphaInd/26)], 'LineWidth', 2)
%     end
    view([0 -90])
    grid
    hold off
    
    axes('Position',[.42 .587 .34 .34])
    hold on
    ylim([-.1 .1])
    xlim([4 9])
    zlim([-2.5 2.5])
    for alphaInd = 1:26
        plot3(Polarization.Alpha(alphaInd) * ones(201, 1), Polarization.Kappa, Polarization.Polarization(alphaInd, :, :), 'Color', [alphaInd/26 0.5 (1 - alphaInd/26)], 'LineWidth', 2)
    end
    view([130 30])
    hold off
% Colors:
    %colormap turbo;
    %c                   = colorbar('eastoutside'); %, 'Direction','reverse');
    %c.Label.String      = '$\left| \Psi  \right|^2$';
    %c.Label.Interpreter = 'latex';
    %c.Label.FontSize    = FontSize + 5; 
% Additional Labels:
    
    %text( 0, 45, 0, '$\epsilon < 0$','HorizontalAlignment', 'center', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize);
% PaperSize:
    set(gcf,'paperunits','in');
    set(gcf,'papersize',[Position(3) + 1, Position(4) + 1]);
% Saving:
    fname   = sprintf('Fig_Polarization_3D_theor.pdf');
    hfig    = gcf;
    print(hfig,'-bestfit','-dpdf', '-r960', fname);