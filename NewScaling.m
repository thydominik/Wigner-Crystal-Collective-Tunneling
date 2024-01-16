%% Figure 5: Spectral gaps for 1,3,5 and 7 particles EXPERIMENTAL
 
clear all
close all
clc

FontSize = 18;
Position = [1 1 7 4];

Figure = figure('color', 'white', 'units', 'inches', 'Position', Position);

hold on
% Fontsize:
    set(gca, 'FontSize', FontSize)

% Axis limits:
    %xlim([0 1100])
    %xticks([-200 200 600 1000])
    ylim([10^-2 50/8])
    %yticks([0 20 40 60])
    %zlim([])
    %axis square
    set(gca, 'Yscale', 'log')
    %set(gca, 'Xscale', 'log')
    %grid
    box on
    scalingVec(1, :) = [7.182 0.875 0.917 0.482]; % pascu's 1st
    scalingVec(2, :) = [13.05 1.327 0.682 0.51];    % Pascu's 2nd
    scalingVec(3, :) = [12.8 9.2 9.067 8.02];   % Original from Shahal
    scalingVec(4, :) = [9 4 3.5 3.5];   % Eyeballed

    Exponent    = 1;
    Factor      = 1;

    
% Axis Labels:
    xlabel('$\tilde \alpha$', 'Interpreter', 'latex', 'FontSize', FontSize+5)
    ylabel('$\Delta^{(N)}(\tilde \alpha)/\Delta^{(N)}(\tilde \alpha = 0)$', 'FontWeight', 'Bold', 'Interpreter', 'latex', 'FontWeight', 'Bold', 'FontSize', FontSize+4)
    %zlabel('$P(\alpha, \epsilon)$', 'Interpreter', 'latex', 'FontSize', FontSize)
% Data: 
    ED1 = load('EDSplitting_1_particle_Nx_5000.mat'); ED1 = ED1.data;
    ED3 = load('EDSplitting_3_particles_restricted_Nx_90_beta_1e-05.mat'); ED3 = ED3.data;
    ED5 = load('ED5P.mat'); ED5 = ED5.data;
    ED7 = load('ED7P.mat'); ED7 = ED7.ED7;
    Exp1 = load('P1.tsv');
    Exp3 = load('P3.tsv');
    Exp5 = load('P5.tsv');
    Exp7 = load('P7.tsv');
% Plotting:
    xshift = 0;
    xscale = 1;
    yscale = 1/1;
    scaleInd = 3;
    fp = 37;
    fx = 0.50 * fp;
    fy = 30;
    Norm = 8;
    %xline(0)

%     plot(Exp1(:, 1)-109, Exp1(:, 2) / Norm, 'o', 'MarkerSize', 8, 'LineWidth', 0.5, MarkerEdgeColor='red', MarkerFaceColor='white');
%     plot((ED1(:, 1) - 2.2) * fx*0.8, ED1(:, 2) * fy / Norm, 's', 'MarkerSize', 8, 'LineWidth', 0.5, MarkerEdgeColor='red', MarkerFaceColor='red')
% 
%     plot((Exp3(:, 1)-(280))*0.25 , Exp3(:, 2) / Norm, 'o', 'MarkerSize', 8, 'LineWidth', 0.5, MarkerEdgeColor='blue', MarkerFaceColor='white')
%     plot((ED3(:, 1) - 6.9) * fx*1.0, ED3(:, 2) * fy / Norm, 's', 'MarkerSize', 8, 'LineWidth', 0.5, MarkerEdgeColor='blue', MarkerFaceColor='blue')
% 
%     plot((Exp5(:, 1)-(560)) * 0.15, Exp5(:, 2) / Norm, 'o', 'MarkerSize', 8, 'LineWidth', 0.5, MarkerEdgeColor='black', MarkerFaceColor='white')
%     plot((ED5(:, 1) -10.4) * fx*1.15, ED5(:, 2) * fy / Norm, 's', 'MarkerSize', 8, 'LineWidth', 0.5, MarkerEdgeColor='black', MarkerFaceColor='black')
% 
%     plot((Exp7(:, 1)-(850))*0.13, Exp7(:, 2) / Norm, 'o', 'MarkerSize', 8, 'LineWidth', 0.5, MarkerEdgeColor='green', MarkerFaceColor='white')
%     plot((ED7(1:end-2, 1) - 13.2) * fx * 1.2, ED7(1:end-2, 2) * fy / Norm, 's', 'MarkerSize', 8, 'LineWidth', 0.5, MarkerEdgeColor='green', MarkerFaceColor='green')
%   



    plot((Exp1(:, 1)-110)/fx, Exp1(:, 2) / Norm, 'o', 'MarkerSize', 8, 'LineWidth', 0.5, MarkerEdgeColor='red', MarkerFaceColor='white');  
    plot((Exp3(:, 1)-(285))*0.25/fx , Exp3(:, 2) / Norm, 'o', 'MarkerSize', 8, 'LineWidth', 0.5, MarkerEdgeColor='blue', MarkerFaceColor='white')
    plot((Exp5(:, 1)-(565)) * 0.15/fx, Exp5(:, 2) / Norm, 'o', 'MarkerSize', 8, 'LineWidth', 0.5, MarkerEdgeColor='black', MarkerFaceColor='white')
    plot((Exp7(:, 1)-(855))*0.13/fx, Exp7(:, 2) / Norm, 'o', 'MarkerSize', 8, 'LineWidth', 0.5, MarkerEdgeColor='green', MarkerFaceColor='white')
    
    plot((ED1(:, 1) - 2.2) *0.8, ED1(:, 2) * fy / Norm, 'o', 'MarkerSize', 5, 'LineWidth', 0.5, MarkerEdgeColor='red', MarkerFaceColor='red')
     plot((ED3(:, 1) - 6.9) *1.0, ED3(:, 2) * fy / Norm, 'o', 'MarkerSize', 5, 'LineWidth', 0.5, MarkerEdgeColor='blue', MarkerFaceColor='blue')
    plot((ED5(:, 1) -10.4) *1.15, ED5(:, 2) * fy / Norm, 'o', 'MarkerSize', 5, 'LineWidth', 0.5, MarkerEdgeColor='black', MarkerFaceColor='black')
    plot((ED7(1:end-2, 1) - 13.2)  * 1.2, ED7(1:end-2, 2) * fy / Norm, 'o', 'MarkerSize', 5, 'LineWidth', 0.5, MarkerEdgeColor='green', MarkerFaceColor='green')
  

% Colors:
    %colormap turbo;
    %c                   = colorbar('eastoutside'); %, 'Direction','reverse');
    %c.Label.String      = '$\left| \Psi  \right|^2$';
    %c.Label.Interpreter = 'latex';
    %c.Label.FontSize    = FontSize + 5; 
% Additional Labels:
    text( -3.75, 1.45*10^-2, 10, '(c)', 'HorizontalAlignment', 'center', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize+4);
    %text( 80, 0.5*10^-1, 0, '$1e^-$', 'HorizontalAlignment', 'center', 'FontWeight', 'Bold', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize);
    %text( 300, 0.5*10^-1, 0, '$3e^-$', 'HorizontalAlignment', 'center', 'FontWeight', 'Bold', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize);
    %text( 600, 0.5*10^-1, 0, '$5e^-$', 'HorizontalAlignment', 'center', 'FontWeight', 'Bold', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize);
    %text( 900, 0.5*10^-1, 0, '$7e^-$', 'HorizontalAlignment', 'center', 'FontWeight', 'Bold', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize);

    qw{1} = plot(nan, 'o', 'MarkerSize', 8, 'LineWidth', 0.5, MarkerEdgeColor='black', MarkerFaceColor='white')
    qw{2} = plot(nan, 'o', 'MarkerSize', 5, 'LineWidth', 0.5, MarkerEdgeColor='black', MarkerFaceColor='black')
    l=legend([qw{:}], {'EXP', 'ED'}, 'location', 'best');
    l.Position=[0.2903 0.3640 0.1429 0.1545];
% PaperSize:
    set(gcf,'paperunits','in');
    set(gcf,'papersize',[Position(3) + 0.2, Position(4) + 0.2]);
% Saving:
    fname   = sprintf('Fig_spectral_gap_exp_log_lin_shahal_scaling_log_v2.pdf');
    hfig    = gcf;
    print(hfig,'-bestfit','-dpdf', '-r960', fname);