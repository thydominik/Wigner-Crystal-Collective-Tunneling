% 
clear all
close all
clc

FontSize = 20;
Position = [1 1 7 6];

Figure = figure('color', 'white', 'units', 'inches', 'Position', Position);

hold on
% Fontsize:
    set(gca, 'FontSize', FontSize)

% Axis limits:
    %xlim([0 1100])
    %xticks([-200 200 600 1000])

    %axis square

    %set(gca, 'Xscale', 'log')
    %grid
    box on
    scalingVec(1, :) = [7.182 0.875 0.917 0.482]; % pascu's 1st
    scalingVec(2, :) = [13.05 1.327 0.682 0.51];    % Pascu's 2nd
    scalingVec(3, :) = [12.8 9.2 9.067 8.02];   % Original from Shahal
    scalingVec(4, :) = [9 4 3.5 3.5];   % Eyeballed

    Exponent    = 1; %1/2.1; % 1.9 - 2.2 majdnem ugyanúgy szemre jót kapunk.
    Factor      = 1;

    %xline([114 368 711 1031])
    %yline([7.182 0.875 0.917 0.482])
    %xline([104 360 720 1020])
    %yline([13.05 1.327 0.682 0.51])
    %xline([104 280 560 850])
    %yline([13.05 9.2 9.067 8.2])


% Data: 
    Exp1 = load('P1.tsv');
    Exp3 = load('P3.tsv');
    Exp5 = load('P5.tsv');
    Exp7 = load('P7.tsv');

    DMRG1   = load('DMRG_Ne_1_gauss.dat');
    DMRG3   = load('DMRG3.mat'); DMRG3 = DMRG3.DMRG3
    DMRG5   = load('DMRG5.mat'); DMRG5 = DMRG5.DMRG5
    DMRG7   = load('DMRG7.mat'); DMRG7 = DMRG7.DMRG7

    ED1 = load('EDSplitting_1_particle_Nx_5000.mat'); ED1 = ED1.data;
    ED3 = load('EDSplitting_3_particles_restricted_Nx_90_beta_1e-05.mat'); ED3 = ED3.data;
    ED5 = load('ED5P.mat'); ED5 = ED5.data;
    ED7 = load('ED7P.mat'); ED7 = ED7.ED7;

    IT1 = load('Standard1particleSplitting.mat'); IT1 = IT1.OneParticleInstanton;
    IT3 = load('Standar3particledSplitting.mat'); IT3 = IT3.SPLITTINGS;
    IT5 = load('Standard5particleSplitting.mat'); IT5 = IT5.SPLITTINGS;
    IT7e = load('Standard7particleSplitting.mat');
    IT7 = IT7e.SPLITTINGS; %IT7e.Instanton7;


% Plotting:
subplot(2,1,1)
hold on
    ylim([10^-3 1*5.57])
    %yticks([10^-3 10 ^-2 10^-1 10^0])
    xlim([0 16])
% Axis Labels:
    xlabel('$\alpha$', 'Interpreter', 'latex', 'FontSize', FontSize)
    ylabel('$\Delta [K]$', 'FontWeight', 'Bold', 'Interpreter', 'latex', 'FontWeight', 'Bold', 'FontSize', FontSize)
    %zlabel('$P(\alpha, \epsilon)$', 'Interpreter', 'latex', 'FontSize', FontSize)
    E = 5.57; % Converted 0.48 meV (energy units) to Kelvin

    scaleInd = 3;
    fy = 30;
    plot((((Exp1(1:end-3, 1) - 109)*1)/14.8) + 2.2        , E * Exp1(1:end-3, 2)/fy, 'x', 'Color', [0.9,0.7,0.0], 'MarkerSize', 8, 'LineWidth', 1.5)
    plot((((Exp3(:, 1) - 280)*0.25)/18.5) + 6.8     , E * Exp3(:, 2)/fy, 'x', 'Color', [0.9,0.7,0.0], 'MarkerSize', 8, 'LineWidth', 1.5)
    plot((((Exp5(:, 1) - 560)*0.15)/21.28) + 10.25  , E * Exp5(:, 2)/fy, 'x', 'Color', [0.9,0.7,0.0], 'MarkerSize', 8, 'LineWidth', 1.5)
    plot((((Exp7(:, 1) - 850)*0.13)/22.2) + 13.05   , E * Exp7(:, 2)/fy, 'x', 'Color', [0.9,0.7,0.0], 'MarkerSize', 8, 'LineWidth', 1.5)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot(ED1(:, 1)              , E * ED1(:, 2)                 , 'rs', 'MarkerSize', 10, 'LineWidth', 1.5)
    plot(-DMRG3(10:end-9, 3)    , E * abs(DMRG3(10:end-9, 6))   , 'rs', 'MarkerSize', 10, 'LineWidth', 1.5)
    plot(-DMRG5(1:end-17, 3)+0.1, E * DMRG5(1:end-17, 6)        , 'rs', 'MarkerSize', 10, 'LineWidth', 1.5)
    plot(DMRG7(1, 1:end-2)+0.2  , E * DMRG7(2, 1:end-2)         , 'rs', 'MarkerSize', 10, 'LineWidth', 1.5)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    plot(ED1(:, 1)      , E * ED1(:, 2)         , 'b.-', 'MarkerSize', 15)
    plot(ED3(:, 1)      , E * ED3(:, 2)         , 'b.-', 'MarkerSize', 15)
    plot(ED5(:, 1)      , E * ED5(:, 2)         , 'b.-', 'MarkerSize', 15)
    plot(ED7(1:end-2, 1), E * ED7(1:end-2, 2)   , 'b.-', 'MarkerSize', 15)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot(IT1(3:end, 1), E * IT1(3:end, 2), 'ko--', 'MarkerFaceColor', 'k')
    plot(IT3(3:end, 1), E * IT3(3:end, 3), 'ko--', 'MarkerFaceColor', 'k')
    plot(IT5(3:end, 1), E * IT5(3:end, 3), 'ko--', 'MarkerFaceColor', 'k')
    plot(IT7(2:end, 1) + 0.2, E * IT7(2:end, 2), 'ko--', 'MarkerFaceColor', 'k')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Colors:
    %colormap turbo;
    %c                   = colorbar('eastoutside'); %, 'Direction','reverse');
    %c.Label.String      = '$\left| \Psi  \right|^2$';
    %c.Label.Interpreter = 'latex';
    %c.Label.FontSize    = FontSize + 5; 
% Additional Labels:
    text( 0.6, 5 * 10^-1, 0, '(a)', 'HorizontalAlignment', 'center', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize);
    text( 4, 2, 0, 'N = 1', 'HorizontalAlignment', 'center', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize-5);
    text( 8.2, 2, 0, 'N = 3', 'HorizontalAlignment', 'center', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize-5);
    text( 11.5, 2, 0, 'N = 5', 'HorizontalAlignment', 'center', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize-5);
    text( 14.2, 2, 0, 'N = 7', 'HorizontalAlignment', 'center', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize-5);
    %text( 80, 0.5*10^-1, 0, '$1e^-$', 'HorizontalAlignment', 'center', 'FontWeight', 'Bold', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize);
    %text( 300, 0.5*10^-1, 0, '$3e^-$', 'HorizontalAlignment', 'center', 'FontWeight', 'Bold', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize);
    %text( 600, 0.5*10^-1, 0, '$5e^-$', 'HorizontalAlignment', 'center', 'FontWeight', 'Bold', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize);
    %text( 900, 0.5*10^-1, 0, '$7e^-$', 'HorizontalAlignment', 'center', 'FontWeight', 'Bold', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize);

    qw{1} = plot(nan, 'x', 'Color', [0.9,0.7,0.0], 'MarkerSize', 8, 'LineWidth', 1.5)
    qw{4} = plot(nan, 'ko--', 'MarkerFaceColor', 'k');
    qw{2} = plot(nan, 'b.-', 'MarkerSize', 15);
    qw{3} = plot(nan, 'rs', 'MarkerSize', 10, 'LineWidth', 1.5)
    legend([qw{:}], {'Exp.', 'ED.', 'DMRG', 'IT.'}, 'FontSize', 16, 'position', [0.91 0.85 0 0])
% PaperSize:
    set(gcf,'paperunits','in');
    set(gcf,'papersize',[Position(3) + 1, Position(4) + 1]);
    
    
    box
    hold off
    subplot(2,1,2)
   
    hold on
    box
        set(gca, 'Yscale', 'log')
        ylim([10^-3 1*5.57])
    yticks([10^-3 10 ^-2 10^-1 10^0])
    xlim([0 16])
    % Axis Labels:
    xlabel('$\alpha$', 'Interpreter', 'latex', 'FontSize', FontSize)
    ylabel('$\Delta [K]$', 'FontWeight', 'Bold', 'Interpreter', 'latex', 'FontWeight', 'Bold', 'FontSize', FontSize)
    %zlabel('$P(\alpha, \epsilon)$', 'Interpreter', 'latex', 'FontSize', FontSize)
    E = 5.57; % Converted 0.48 meV (energy units) to Kelvin

    scaleInd = 3;
    fy = 30;
    plot((((Exp1(1:end-3, 1) - 109)*1)/14.8) + 2.2        , E * Exp1(1:end-3, 2)/fy, 'x', 'Color', [0.9,0.7,0.0], 'MarkerSize', 8, 'LineWidth', 1.5)
    plot((((Exp3(:, 1) - 280)*0.25)/18.5) + 6.8     , E * Exp3(:, 2)/fy, 'x', 'Color', [0.9,0.7,0.0], 'MarkerSize', 8, 'LineWidth', 1.5)
    plot((((Exp5(:, 1) - 560)*0.15)/21.28) + 10.25  , E * Exp5(:, 2)/fy, 'x', 'Color', [0.9,0.7,0.0], 'MarkerSize', 8, 'LineWidth', 1.5)
    plot((((Exp7(:, 1) - 850)*0.13)/22.2) + 13.05   , E * Exp7(:, 2)/fy, 'x', 'Color', [0.9,0.7,0.0], 'MarkerSize', 8, 'LineWidth', 1.5)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot(ED1(:, 1)              , E * ED1(:, 2)                 , 'rs', 'MarkerSize', 10, 'LineWidth', 1.5)
    plot(-DMRG3(10:end-9, 3)    , E * abs(DMRG3(10:end-9, 6))   , 'rs', 'MarkerSize', 10, 'LineWidth', 1.5)
    plot(-DMRG5(1:end-17, 3)+0.1, E * DMRG5(1:end-17, 6)        , 'rs', 'MarkerSize', 10, 'LineWidth', 1.5)
    plot(DMRG7(1, 1:end-2)+0.2  , E * DMRG7(2, 1:end-2)         , 'rs', 'MarkerSize', 10, 'LineWidth', 1.5)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    plot(ED1(3:end, 1)      , E * ED1(3:end, 2)         , 'b.-', 'MarkerSize', 15)
    plot(ED3(:, 1)      , E * ED3(:, 2)         , 'b.-', 'MarkerSize', 15)
    plot(ED5(:, 1)      , E * ED5(:, 2)         , 'b.-', 'MarkerSize', 15)
    plot(ED7(1:end-2, 1), E * ED7(1:end-2, 2)   , 'b.-', 'MarkerSize', 15)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot(IT1(3:end, 1), E * IT1(3:end, 2), 'ko--', 'MarkerFaceColor', 'k')
    plot(IT3(3:end, 1), E * IT3(3:end, 3), 'ko--', 'MarkerFaceColor', 'k')
    plot(IT5(3:end, 1), E * IT5(3:end, 3), 'ko--', 'MarkerFaceColor', 'k')
    plot(IT7(2:end, 1) + 0.2, E * IT7(2:end, 2), 'ko--', 'MarkerFaceColor', 'k')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Colors:
    %colormap turbo;
    %c                   = colorbar('eastoutside'); %, 'Direction','reverse');
    %c.Label.String      = '$\left| \Psi  \right|^2$';
    %c.Label.Interpreter = 'latex';
    %c.Label.FontSize    = FontSize + 5; 
% Additional Labels:
    text(  0.6, .22 * 10^-2, 0, '(b)', 'HorizontalAlignment', 'center', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize);
    %text( 4, 2, 0, 'N = 1', 'HorizontalAlignment', 'center', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize-5);
    %text( 8.2, 2, 0, 'N = 3', 'HorizontalAlignment', 'center', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize-5);
    %text( 11.5, 2, 0, 'N = 5', 'HorizontalAlignment', 'center', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize-5);
    %text( 14.2, 2, 0, 'N = 7', 'HorizontalAlignment', 'center', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize-5);
    %text( 80, 0.5*10^-1, 0, '$1e^-$', 'HorizontalAlignment', 'center', 'FontWeight', 'Bold', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize);
    %text( 300, 0.5*10^-1, 0, '$3e^-$', 'HorizontalAlignment', 'center', 'FontWeight', 'Bold', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize);
    %text( 600, 0.5*10^-1, 0, '$5e^-$', 'HorizontalAlignment', 'center', 'FontWeight', 'Bold', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize);
    %text( 900, 0.5*10^-1, 0, '$7e^-$', 'HorizontalAlignment', 'center', 'FontWeight', 'Bold', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize);

    qw{1} = plot(nan, 'x', 'Color', [0.9,0.7,0.0], 'MarkerSize', 8, 'LineWidth', 1.5)
    qw{4} = plot(nan, 'ko--', 'MarkerFaceColor', 'k');
    qw{2} = plot(nan, 'b.-', 'MarkerSize', 15);
    qw{3} = plot(nan, 'rs', 'MarkerSize', 10, 'LineWidth', 1.5)
    %legend([qw{:}], {'Exp.', 'ED.', 'DMRG', 'IT.'}, 'FontSize', 16, 'position', [0.91 0.8 0 0])
% PaperSize:
    set(gcf,'paperunits','in');
    set(gcf,'papersize',[Position(3) + 1, Position(4) + 1]);
    hold off
% Saving:
    fname   = sprintf('Fig_spectral_gap_exp_fitted_2.pdf');
    hfig    = gcf;
    print(hfig,'-bestfit','-dpdf', '-r960', fname);