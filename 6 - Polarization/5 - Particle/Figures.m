clc
clear all
close all
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
func = WF1 + WF2 + WF3 + WF4 + WF5;
Int = 0;
for i = 2:length(func)
    Int = Int + 160 * abs(x(i) - x(i-1)) * 0.5 * (func(i) + func(i-1));
end
func = func/Int * 5
plot(x * 160/1000, func, 'r-', 'LineWidth', 6)
x = linspace(-7, 7, 1000);
LinePlot = plot(x*160 / 1000, 3/160 * 0.1 * 10^-2 * (0.25 * (x.^2 - WF.alpha).^2 + WF.kappa * x), 'k-', 'LineWidth', 5);
LinePlot.Color(4) = 0.5;
set(gca, 'FontSize', FontSize)
xlim([-7*160/1000 7*160/1000])
xticks([-1 0 1])
ylim([-0.001 0.01])
yticks([0 0.01])
xlabel('$z [\mu{\rm{m}}]$', 'Interpreter', 'latex', 'FontSize', FontSize)
label = ylabel('$\rho(z) [\mu{\rm{m}}^{-1}] $', 'FontWeight', 'Bold', 'Interpreter', 'latex', 'FontWeight', 'Bold', 'FontSize', FontSize + 5, 'Color','r')
%yline(0, 'k-', 'LineWidth', 2.5, 'Alpha', 1);
text( 0, 1.4/160, 0, '$\epsilon = 0$','HorizontalAlignment', 'center', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize + 15);
set(gcf,'paperunits','in');
set(gcf,'papersize',[Position(3) + 1, Position(4) + 1]);
box
hold off
fname   = sprintf('Fig_WaveFunction_Potential_eps1_v2.pdf');
hfig    = gcf;
print(hfig,'-bestfit','-dpdf', '-r960', fname);
%%

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
func = WF1 + WF2 + WF3 + WF4 + WF5;
Int = 0;
for i = 2:length(func)
    Int = Int + 160 * abs(x(i) - x(i-1)) * 0.5 * (func(i) + func(i-1));
end
func = func/Int * 5
plot(x * 160/1000, func, 'r-', 'LineWidth', 6)
x = linspace(-7, 7, 1000);
LinePlot = plot(x*160/1000, 3/160 * 0.1 * 10^-2 * (0.25 * (x.^2 - WF.alpha).^2 + WF.kappa * x), 'k-', 'LineWidth', 5);
LinePlot.Color(4) = 0.5;
set(gca, 'FontSize', FontSize)
xlim([-7*160/1000 7*160/1000])
xticks([-1 0 1])
ylim([-0.001 0.01])
yticks([0 0.01])
xlabel('$z [\mu{\rm{m}}]$', 'Interpreter', 'latex', 'FontSize', FontSize)
label = ylabel('$\rho(z) [\mu{\rm{m}}^{-1}] $', 'FontWeight', 'Bold', 'Interpreter', 'latex', 'FontWeight', 'Bold', 'FontSize', FontSize + 5, 'Color','r')
%yline(0, 'k-', 'LineWidth', 2.5, 'Alpha', 1);
text( 0, 1.37/160, 0, '$\epsilon = -10^{-6}$','HorizontalAlignment', 'center', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize + 15);
set(gcf,'paperunits','in');
set(gcf,'papersize',[Position(3) + 1, Position(4) + 1]);
box
hold off
fname   = sprintf('Fig_WaveFunction_Potential_eps2_v2.pdf');
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
func = WF1 + WF2 + WF3 + WF4 + WF5;
Int = 0;
for i = 2:length(func)
    Int = Int + 160 * abs(x(i) - x(i-1)) * 0.5 * (func(i) + func(i-1));
end
func = func/Int * 5
plot(x * 160/1000, func, 'r-', 'LineWidth', 6)
x = linspace(-7, 7, 1000);
LinePlot = plot(x*160/1000, 3/160 * 0.1 * 10^-2 * (0.25 * (x.^2 - WF.alpha).^2 + WF.kappa * x), 'k-', 'LineWidth', 5);
LinePlot.Color(4) = 0.5;
set(gca, 'FontSize', FontSize)
xlim([-7*160/1000 7*160/1000])
xticks([-1 0 1])
ylim([-0.001 0.01])
yticks([0 0.01])
xlabel('$z [\mu{\rm{m}}]$', 'Interpreter', 'latex', 'FontSize', FontSize)
label = ylabel('$\rho(z) [\mu{\rm{m}}^{-1}] $', 'FontWeight', 'Bold', 'Interpreter', 'latex', 'FontWeight', 'Bold', 'FontSize', FontSize + 5, 'Color','r')
%yline(0, 'k-', 'LineWidth', 2.5, 'Alpha', 1);
text( 0, 1.37/160, 0, '$\epsilon = -2\cdot 10^{-6}$','HorizontalAlignment', 'center', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize + 15);
set(gcf,'paperunits','in');
set(gcf,'papersize',[Position(3) + 1, Position(4) + 1]);
box
hold off
fname   = sprintf('Fig_WaveFunction_Potential_eps3_v2.pdf');
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
func = WF1 + WF2 + WF3 + WF4 + WF5;
Int = 0;
for i = 2:length(func)
    Int = Int + 160 * abs(x(i) - x(i-1)) * 0.5 * (func(i) + func(i-1));
end
func = func/Int * 5
plot(x * 160/1000, func, 'r-', 'LineWidth', 6)
x = linspace(-7, 7, 1000);
LinePlot = plot(x*160/1000, 3/160 * 0.1 * 10^-2 * (0.25 * (x.^2 - WF.alpha).^2 + WF.kappa * x), 'k-', 'LineWidth', 5);
LinePlot.Color(4) = 0.5;
set(gca, 'FontSize', FontSize)
xlim([-7*160/1000 7*160/1000])
xticks([-1 0 1])
ylim([-0.001 0.01])
yticks([0 0.01])
xlabel('$z [\mu{\rm{m}}]$', 'Interpreter', 'latex', 'FontSize', FontSize)
label = ylabel('$\rho(z) [\mu{\rm{m}}^{-1}] $', 'FontWeight', 'Bold', 'Interpreter', 'latex', 'FontWeight', 'Bold', 'FontSize', FontSize + 5, 'Color','r')
%yline(0, 'k-', 'LineWidth', 2.5, 'Alpha', 1);
text( 0, 1.37/160, 0, '$\epsilon = -5\cdot 10^{-6}$','HorizontalAlignment', 'center', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize + 15);
set(gcf,'paperunits','in');
set(gcf,'papersize',[Position(3) + 1, Position(4) + 1]);
box
hold off
fname   = sprintf('Fig_WaveFunction_Potential_eps4_v2.pdf');
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
label = ylabel('$\rho(z) [{\rm{nm}}^{-1}] $', 'FontWeight', 'Bold', 'Interpreter', 'latex', 'FontWeight', 'Bold', 'FontSize', FontSize + 5, 'Color','r')
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
WF = load("WaveFunction5e-05.mat"); WF = WF.WF;
x = linspace(-7, 7, 100);
WF1 = interp1(WF.x1, WF.WF1, x, 'cubic'); WF1 ( isnan(WF1)) = 0; WF1 = WF1 / norm(WF1);
WF2 = interp1(WF.x2, WF.WF2, x, 'cubic'); WF2 ( isnan(WF2)) = 0;WF2 = WF2 / norm(WF2);
WF3 = interp1(WF.x3, WF.WF3, x, 'cubic'); WF3 ( isnan(WF3)) = 0;WF3 = WF3 / norm(WF3);
WF4 = interp1(WF.x4, WF.WF4, x, 'cubic'); WF4 ( isnan(WF4)) = 0;WF4 = WF4 / norm(WF4);
WF5 = interp1(WF.x5, WF.WF5, x, 'cubic'); WF5 ( isnan(WF5)) = 0;WF5 = WF5 / norm(WF5);
f =  -(WF1 + WF2 + WF3 + WF4 + WF5)/norm((WF1 + WF2 + WF3 + WF4 +WF5));
a = 0.2;
b = 0.03
c = 1
C = conv(f - flip(f), (normpdf([-a:b:a],0,c)));
plot(linspace(-7, 7, length(C)) * 160/1000,  - C , 'r-', 'LineWidth', 6)
%plot(x * 160,  -(f - flip(f)), 'b-', 'LineWidth', 3)
%plot(x*160, 10*[WF1;WF2; WF3; WF4; WF5], 'r.')
%plot(x*160, 10*[flip(WF1);flip(WF2); flip(WF3); flip(WF4); flip(WF5)], 'bo')
%plot([-a:b:a] * 160, normpdf([-a:b:a],0,c))
x = linspace(-7, 7, 1000);
%LinePlot = plot(x*160, 0.1 * 10^-2 * (0.25 * (x.^2 - WF.alpha).^2 + WF.kappa * x), 'k-', 'LineWidth', 5);
%LinePlot.Color(4) = 0.5;
set(gca, 'FontSize', FontSize)
%xlim([-7*160 7*160])
%xticks([-500 0 500])
%ylim([-0.8 0.8])
yticks([-1 1])
xlabel('$z [\mu{\rm{m}}]$', 'Interpreter', 'latex', 'FontSize', FontSize)
%label = ylabel('$\rho(z) \left[ \frac{1}{nm}\right] $', 'FontWeight', 'Bold', 'Interpreter', 'latex', 'FontWeight', 'Bold', 'FontSize', FontSize + 5, 'Color','r')
x = linspace(-7, 7, 1000);
mult = 100;
yline(0, 'k-', 'LineWidth', 2.5, 'Alpha', 1);
trshld = 8*10^-4;
shift1 = -9.5;
text( 0.5, 0.3, 0, '$N = 5$','HorizontalAlignment', 'center', 'Color', 'black', 'interpreter', 'Latex', 'FontSize', FontSize + 15);
set(gcf,'paperunits','in');
set(gcf,'papersize',[Position(3) + 1, Position(4) + 1]);
box
hold off
fname   = sprintf('Fig_WaveFunction_difference_N_5.pdf');
hfig    = gcf;
print(hfig,'-bestfit','-dpdf', '-r960', fname);

%%
clf
ax = tight_subplot(2,2,0.01,0,.1);hold on
set(ax,'color','none')
set(ax(3:4),'xcolor','none')
set(ax([1 3]),'ycolor','none')
set([ax.XRuler],'tickdirection','both')
set([ax.YRuler],'tickdirection','both')
for i=1:4
    axes(ax(i));hold on
    plot(rand(10,10))
end
linkaxes(ax,'xy')



%%

function [ha, pos] = tight_subplot(Nh, Nw, gap, marg_h, marg_w)
% tight_subplot creates "subplot" axes with adjustable gaps and margins
%
% [ha, pos] = tight_subplot(Nh, Nw, gap, marg_h, marg_w)
%
%   in:  Nh      number of axes in hight (vertical direction)
%        Nw      number of axes in width (horizontaldirection)
%        gap     gaps between the axes in normalized units (0...1)
%                   or [gap_h gap_w] for different gaps in height and width 
%        marg_h  margins in height in normalized units (0...1)
%                   or [lower upper] for different lower and upper margins 
%        marg_w  margins in width in normalized units (0...1)
%                   or [left right] for different left and right margins 
%
%  out:  ha     array of handles of the axes objects
%                   starting from upper left corner, going row-wise as in
%                   subplot
%        pos    positions of the axes objects
%
%  Example: ha = tight_subplot(3,2,[.01 .03],[.1 .01],[.01 .01])
%           for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end
%           set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')
% Pekka Kumpulainen 21.5.2012   @tut.fi
% Tampere University of Technology / Automation Science and Engineering
if nargin<3; gap = .02; end
if nargin<4 || isempty(marg_h); marg_h = .05; end
if nargin<5; marg_w = .05; end
if numel(gap)==1; 
    gap = [gap gap];
end
if numel(marg_w)==1; 
    marg_w = [marg_w marg_w];
end
if numel(marg_h)==1; 
    marg_h = [marg_h marg_h];
end
axh = (1-sum(marg_h)-(Nh-1)*gap(1))/Nh; 
axw = (1-sum(marg_w)-(Nw-1)*gap(2))/Nw;
py = 1-marg_h(2)-axh; 
% ha = zeros(Nh*Nw,1);
ii = 0;
for ih = 1:Nh
    px = marg_w(1);
    
    for ix = 1:Nw
        ii = ii+1;
        ha(ii) = axes('Units','normalized', ...
            'Position',[px py axw axh], ...
            'XTickLabel','', ...
            'YTickLabel','');
        px = px+axw+gap(2);
    end
    py = py-axh-gap(1);
end
if nargout > 1
    pos = get(ha,'Position');
end
ha = ha(:);
end