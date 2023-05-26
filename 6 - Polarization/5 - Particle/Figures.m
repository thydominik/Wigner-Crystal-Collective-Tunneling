clc
clear all

FontSize = 30;
Position = [1 1 6 5];

Figure = figure('color', 'white', 'units', 'inches', 'Position', Position);
hold on
WF = load("WaveFunction0.mat"); WF = WF.WF;
x = linspace(-7, 7, 1000);
WF1 = interp1(WF.x1, WF.WF1, x, 'cubic'); WF1 ( isnan(WF1)) = 0; WF1 = WF1 / norm(WF1);
WF2 = interp1(WF.x2, WF.WF2, x, 'cubic'); WF2 ( isnan(WF2)) = 0;WF2 = WF2 / norm(WF2);
WF3 = interp1(WF.x3, WF.WF3, x, 'cubic'); WF3 ( isnan(WF3)) = 0;WF3 = WF3 / norm(WF3);
WF4 = interp1(WF.x4, WF.WF4, x, 'cubic'); WF4 ( isnan(WF4)) = 0;WF4 = WF4 / norm(WF4);
WF5 = interp1(WF.x5, WF.WF5, x, 'cubic'); WF5 ( isnan(WF5)) = 0;WF5 = WF5 / norm(WF5);
plot(x * 160, 5 * abs(0.5 *(WF1 + flip(WF1) + WF2 + flip(WF2) + flip(WF3) + flip(WF4) + flip(WF5) + WF3 + WF4 + WF5))/norm(0.5 *(WF1 + flip(WF1) + WF2 + flip(WF2) + flip(WF3) + flip(WF4) + flip(WF5) + WF3 + WF4 + WF5)), 'r-', 'LineWidth', 6)
x = linspace(-7, 7, 1000);
LinePlot = plot(x*160, 0.1 * 10^-2 * (0.25 * (x.^2 - WF.alpha).^2 + WF.kappa * x), 'k-', 'LineWidth', 5);
LinePlot.Color(4) = 0.5;
set(gca, 'FontSize', FontSize)
xlim([-7*160 7*160])
xticks([-500 0 500])
ylim([-0.02 0.5])
yticks([0 0.5])
xlabel('$z [{\rm{nm}}]$', 'Interpreter', 'latex', 'FontSize', FontSize)
label = ylabel('$\rho(z) \left[ \frac{1}{nm}\right] $', 'FontWeight', 'Bold', 'Interpreter', 'latex', 'FontWeight', 'Bold', 'FontSize', FontSize + 5, 'Color','r')
x = linspace(-7, 7, 1000);
mult = 100;
yline(0, 'k-', 'LineWidth', 2.5, 'Alpha', 1);
trshld = 8*10^-4;
shift1 = -9.5;
text( 0, 0.45, 0, '$\epsilon = 0$','HorizontalAlignment', 'center', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize + 15);
set(gcf,'paperunits','in');
set(gcf,'papersize',[Position(3) + 1, Position(4) + 1]);
box
hold off
fname   = sprintf('Fig_WaveFunction_Potential_eps1.pdf');
hfig    = gcf;
print(hfig,'-bestfit','-dpdf', '-r960', fname);
%%
close all
FontSize = 30;
Position = [1 1 6 5];

Figure = figure('color', 'white', 'units', 'inches', 'Position', Position);
hold on
WF = load("WaveFunction1e-05.mat"); WF = WF.WF;
x = linspace(-7, 7, 1000);
WF1 = interp1(WF.x1, WF.WF1, x, 'cubic'); WF1 ( isnan(WF1)) = 0; WF1 = WF1 / norm(WF1);
WF2 = interp1(WF.x2, WF.WF2, x, 'cubic'); WF2 ( isnan(WF2)) = 0;WF2 = WF2 / norm(WF2);
WF3 = interp1(WF.x3, WF.WF3, x, 'cubic'); WF3 ( isnan(WF3)) = 0;WF3 = WF3 / norm(WF3);
WF4 = interp1(WF.x4, WF.WF4, x, 'cubic'); WF4 ( isnan(WF4)) = 0;WF4 = WF4 / norm(WF4);
WF5 = interp1(WF.x5, WF.WF5, x, 'cubic'); WF5 ( isnan(WF5)) = 0;WF5 = WF5 / norm(WF5);
plot(x * 160, 5 * abs(smooth(WF1 + WF2 + WF3 + WF4 +WF5))/norm(abs(smooth(WF1 + WF2 + WF3 + WF4 +WF5))), 'r-', 'LineWidth', 6)
x = linspace(-7, 7, 1000);
LinePlot = plot(x*160, 0.1 * 10^-2 * (0.25 * (x.^2 - WF.alpha).^2 + WF.kappa * x), 'k-', 'LineWidth', 5);
LinePlot.Color(4) = 0.5;
set(gca, 'FontSize', FontSize)
xlim([-7*160 7*160])
xticks([-500 0 500])
ylim([-0.02 0.5])
yticks([0 0.5])
xlabel('$z [{\rm{nm}}]$', 'Interpreter', 'latex', 'FontSize', FontSize)
label = ylabel('$\rho(z) \left[ \frac{1}{nm}\right] $', 'FontWeight', 'Bold', 'Interpreter', 'latex', 'FontWeight', 'Bold', 'FontSize', FontSize + 5, 'Color','r')
x = linspace(-7, 7, 1000);
mult = 100;
yline(0, 'k-', 'LineWidth', 2.5, 'Alpha', 1);
trshld = 8*10^-4;
shift1 = -9.5;
text( 10, 0.5, 0, '$\epsilon = 10^{-6}$','HorizontalAlignment', 'center', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize + 15);
set(gcf,'paperunits','in');
set(gcf,'papersize',[Position(3) + 1, Position(4) + 1]);
box
hold off
fname   = sprintf('Fig_WaveFunction_Potential_eps2.pdf');
hfig    = gcf;
print(hfig,'-bestfit','-dpdf', '-r960', fname);

%%
close all
FontSize = 30;
Position = [1 1 6 5];

Figure = figure('color', 'white', 'units', 'inches', 'Position', Position);
hold on
WF = load("WaveFunction2e-05.mat"); WF = WF.WF;
x = linspace(-7, 7, 1000);
WF1 = interp1(WF.x1, WF.WF1, x, 'cubic'); WF1 ( isnan(WF1)) = 0; WF1 = WF1 / norm(WF1);
WF2 = interp1(WF.x2, WF.WF2, x, 'cubic'); WF2 ( isnan(WF2)) = 0;WF2 = WF2 / norm(WF2);
WF3 = interp1(WF.x3, WF.WF3, x, 'cubic'); WF3 ( isnan(WF3)) = 0;WF3 = WF3 / norm(WF3);
WF4 = interp1(WF.x4, WF.WF4, x, 'cubic'); WF4 ( isnan(WF4)) = 0;WF4 = WF4 / norm(WF4);
WF5 = interp1(WF.x5, WF.WF5, x, 'cubic'); WF5 ( isnan(WF5)) = 0;WF5 = WF5 / norm(WF5);
plot(x*160, 5 * smooth(WF1 + WF2 + WF3 + WF4 +WF5)/ norm(smooth(WF1 + WF2 + WF3 + WF4 +WF5)), 'r-', 'LineWidth', 6)
x = linspace(-7, 7, 1000);
LinePlot = plot(x*160, 0.1 * 10^-2 * (0.25 * (x.^2 - WF.alpha).^2 + WF.kappa * x), 'k-', 'LineWidth', 5);
LinePlot.Color(4) = 0.5;
set(gca, 'FontSize', FontSize)
xlim([-7*160 7*160])
xticks([-500 0 500])
ylim([-0.02 0.5])
yticks([0 0.5])
xlabel('$z [{\rm{nm}}]$', 'Interpreter', 'latex', 'FontSize', FontSize)
label = ylabel('$\rho(z) \left[ \frac{1}{nm}\right] $', 'FontWeight', 'Bold', 'Interpreter', 'latex', 'FontWeight', 'Bold', 'FontSize', FontSize + 5, 'Color','r')
x = linspace(-7, 7, 1000);
mult = 100;
yline(0, 'k-', 'LineWidth', 2.5, 'Alpha', 1);
trshld = 8*10^-4;
shift1 = -9.5;
text( 10, 0.42, 0, '$\epsilon = 2\cdot10^{-6}$','HorizontalAlignment', 'center', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize + 15);
set(gcf,'paperunits','in');
set(gcf,'papersize',[Position(3) + 1, Position(4) + 1]);
box
hold off
fname   = sprintf('Fig_WaveFunction_Potential_eps3.pdf');
hfig    = gcf;
print(hfig,'-bestfit','-dpdf', '-r960', fname);


%%

close all
FontSize = 30;
Position = [1 1 6 5];

Figure = figure('color', 'white', 'units', 'inches', 'Position', Position);
hold on
WF = load("WaveFunction5e-05.mat"); WF = WF.WF;
x = linspace(-7, 7, 1000);
WF1 = interp1(WF.x1, WF.WF1, x, 'cubic'); WF1 ( isnan(WF1)) = 0; WF1 = WF1 / norm(WF1);
WF2 = interp1(WF.x2, WF.WF2, x, 'cubic'); WF2 ( isnan(WF2)) = 0;WF2 = WF2 / norm(WF2);
WF3 = interp1(WF.x3, WF.WF3, x, 'cubic'); WF3 ( isnan(WF3)) = 0;WF3 = WF3 / norm(WF3);
WF4 = interp1(WF.x4, WF.WF4, x, 'pchip'); WF4 ( isnan(WF4)) = 0;WF4 = WF4 / norm(WF4);
WF5 = interp1(WF.x5, WF.WF5, x, 'cubic'); WF5 ( isnan(WF5)) = 0;WF5 = WF5 / norm(WF5);
plot(x * 160, 5 * smooth(WF1 + WF2 + WF3 + WF4 +WF5)/norm(smooth(WF1 + WF2 + WF3 + WF4 +WF5)), 'r-', 'LineWidth', 6)
x = linspace(-7, 7, 1000);
LinePlot = plot(x*160, 0.1 * 10^-2 * (0.25 * (x.^2 - WF.alpha).^2 + WF.kappa * x), 'k-', 'LineWidth', 5);
LinePlot.Color(4) = 0.5;
set(gca, 'FontSize', FontSize)
xlim([-7*160 7*160])
xticks([-500 0 500])
ylim([-0.02 0.5])
yticks([0 0.5])
xlabel('$z [{\rm{nm}}]$', 'Interpreter', 'latex', 'FontSize', FontSize)
label = ylabel('$\rho(z) \left[ \frac{1}{nm}\right] $', 'FontWeight', 'Bold', 'Interpreter', 'latex', 'FontWeight', 'Bold', 'FontSize', FontSize + 5, 'Color','r')
x = linspace(-7, 7, 1000);
mult = 100;
yline(0, 'k-', 'LineWidth', 2.5, 'Alpha', 1);
trshld = 8*10^-4;
shift1 = -9.5;
text( 10, 0.42, 0, '$\epsilon = 5\cdot10^{-6}$','HorizontalAlignment', 'center', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize + 15);
set(gcf,'paperunits','in');
set(gcf,'papersize',[Position(3) + 1, Position(4) + 1]);
box
hold off
fname   = sprintf('Fig_WaveFunction_Potential_eps4.pdf');
hfig    = gcf;
print(hfig,'-bestfit','-dpdf', '-r960', fname);
%%

close all
FontSize = 30;
Position = [1 1 6 5];

Figure = figure('color', 'white', 'units', 'inches', 'Position', Position);
hold on
WF = load("WaveFunction0_001.mat"); WF = WF.WF;
x = linspace(-7, 7, 1000);
WF1 = interp1(WF.x1, WF.WF1, x, 'cubic'); WF1 ( isnan(WF1)) = 0; WF1 = WF1 / norm(WF1);
WF2 = interp1(WF.x2, WF.WF2, x, 'cubic'); WF2 ( isnan(WF2)) = 0;WF2 = WF2 / norm(WF2);
WF3 = interp1(WF.x3, WF.WF3, x, 'cubic'); WF3 ( isnan(WF3)) = 0;WF3 = WF3 / norm(WF3);
WF4 = interp1(WF.x4, WF.WF4, x, 'pchip'); WF4 ( isnan(WF4)) = 0;WF4 = WF4 / norm(WF4);
WF5 = interp1(WF.x5, WF.WF5, x, 'cubic'); WF5 ( isnan(WF5)) = 0;WF5 = WF5 / norm(WF5);
plot(x * 160, 5 * smooth(WF1 + WF2 + WF3 + WF4 +WF5)/norm(smooth(WF1 + WF2 + WF3 + WF4 +WF5)), 'r-', 'LineWidth', 6)
x = linspace(-7, 7, 1000);
LinePlot = plot(x*160, 0.1 * 10^-2 * (0.25 * (x.^2 - WF.alpha).^2 + WF.kappa * x), 'k-', 'LineWidth', 5);
LinePlot.Color(4) = 0.5;
set(gca, 'FontSize', FontSize)
xlim([-7*160 7*160])
xticks([-500 0 500])
ylim([-0.02 0.5])
yticks([0 0.5])
xlabel('$z [{\rm{nm}}]$', 'Interpreter', 'latex', 'FontSize', FontSize)
label = ylabel('$\rho(z) \left[ \frac{1}{nm}\right] $', 'FontWeight', 'Bold', 'Interpreter', 'latex', 'FontWeight', 'Bold', 'FontSize', FontSize + 5, 'Color','r')
x = linspace(-7, 7, 1000);
mult = 100;
yline(0, 'k-', 'LineWidth', 2.5, 'Alpha', 1);
trshld = 8*10^-4;
shift1 = -9.5;
text( 10, 0.42, 0, '$\epsilon = 10^{-4}$','HorizontalAlignment', 'center', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize + 15);
set(gcf,'paperunits','in');
set(gcf,'papersize',[Position(3) + 1, Position(4) + 1]);
box
hold off
fname   = sprintf('Fig_WaveFunction_Potential_eps5.pdf');
hfig    = gcf;
print(hfig,'-bestfit','-dpdf', '-r960', fname);

%%

close all
FontSize = 30;
Position = [1 1 6 5];

Figure = figure('color', 'white', 'units', 'inches', 'Position', Position);
hold on
WF = load("WaveFunction0_001.mat"); WF = WF.WF;
x = linspace(-7, 7, 1000);
WF1 = interp1(WF.x1, WF.WF1, x, 'cubic'); WF1 ( isnan(WF1)) = 0; WF1 = WF1 / norm(WF1);
WF2 = interp1(WF.x2, WF.WF2, x, 'cubic'); WF2 ( isnan(WF2)) = 0;WF2 = WF2 / norm(WF2);
WF3 = interp1(WF.x3, WF.WF3, x, 'cubic'); WF3 ( isnan(WF3)) = 0;WF3 = WF3 / norm(WF3);
WF4 = interp1(WF.x4, WF.WF4, x, 'pchip'); WF4 ( isnan(WF4)) = 0;WF4 = WF4 / norm(WF4);
WF5 = interp1(WF.x5, WF.WF5, x, 'cubic'); WF5 ( isnan(WF5)) = 0;WF5 = WF5 / norm(WF5);
f = 5 * (WF1 + WF2 + WF3 + WF4 +WF5)/norm((WF1 + WF2 + WF3 + WF4 +WF5));
C = conv(f - flip(f), normpdf([-2:0.02:2],0,1));
plot(linspace(-7, 7, length(C)), 10^-2 * C , 'r-', 'LineWidth', 6)
plot(x, WF1)
plot(x, WF2)
plot(x, WF3)
plot(x, WF4)
plot(x, WF5)
x = linspace(-7, 7, 1000);
%LinePlot = plot(x*160, 0.1 * 10^-2 * (0.25 * (x.^2 - WF.alpha).^2 + WF.kappa * x), 'k-', 'LineWidth', 5);
%LinePlot.Color(4) = 0.5;
set(gca, 'FontSize', FontSize)
%xlim([-7*160 7*160])
%xticks([-500 0 500])
%ylim([-0.02 0.5])
%yticks([0 0.5])
xlabel('$z [{\rm{nm}}]$', 'Interpreter', 'latex', 'FontSize', FontSize)
label = ylabel('$\rho(z) \left[ \frac{1}{nm}\right] $', 'FontWeight', 'Bold', 'Interpreter', 'latex', 'FontWeight', 'Bold', 'FontSize', FontSize + 5, 'Color','r')
x = linspace(-7, 7, 1000);
mult = 100;
yline(0, 'k-', 'LineWidth', 2.5, 'Alpha', 1);
trshld = 8*10^-4;
shift1 = -9.5;
%text( 10, 0.42, 0, '$\epsilon = 10^{-4}$','HorizontalAlignment', 'center', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize + 15);
set(gcf,'paperunits','in');
set(gcf,'papersize',[Position(3) + 1, Position(4) + 1]);
box
hold off
fname   = sprintf('Fig_WaveFunction_difference_foreps_e-4.pdf');
hfig    = gcf;
print(hfig,'-bestfit','-dpdf', '-r960', fname);