%% General cleaning

close all;
clear all;
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultFigureWindowStyle','normal');



%% Generates the 2d density plot for density average for U=0, ballistic.
close all;
clear all;

mu = 0.0005;
L = 200;
Mmax = 200;
U = 0.00;
dtm = 0.2;

datafile = ['./TEBD_no_leads_L_', num2str(L),'_U_',num2str(U, '%.2f'), '_mu_', num2str(mu, '%.4f'), '_Mmax_', num2str(Mmax) , '/datafile.mat']
load(datafile);

average_occupation_TEBD = (average_occupation_TEBD-1)/(mu/2);
h2=figure('color','white','units','inches','position',[1 1 7 5]);
hold on; 

size_x = size(average_occupation_TEBD, 2);
size_t = size(average_occupation_TEBD, 1);

sites = [-size_x/2:size_x/2-1];
times = [1:size_t]*dtm;
pcolor(sites, times, average_occupation_TEBD);

shading interp;

xlim ([-78,78]);
ylim([0,50]);

fsize = 20;
set(gca, 'Fontsize', fsize);
xlabel(sprintf('$x$'),'Interpreter', 'latex', 'FontSize',fsize+5);
ylabel(sprintf('$t$'),'Interpreter', 'latex', 'FontSize',fsize+5);


x = [0:-0.1:-70];
y = -x;
plot(x,y,'Color','r','LineStyle','--')
x = [0:0.1:70];
y = x;
line(x,y,'Color','r','LineStyle','--')



% find the <n> =0.5 and plot it as lines. 
% 
% n1_half = [];
% n2_half = [];
% 
% for t=5:size_t
%     [x, idx] = unique(average_occupation_TEBD(t, :));
%     
%     y = sites(idx);
% 
%     x_half = interp1(x,y, 0.5);
%     n1_half =[n1_half; x_half, times(t)]; 
% 
%     x_half = interp1(x, y, -0.5);
%     n2_half =[n2_half; x_half, times(t)]; 
% end
% 
% plot(n1_half(1:10:end,1), smooth(n1_half(1:10:end,2)), 'ro' );
% plot(n2_half(1:10:end,1), n2_half(1:10:end,2), 'ro' );
% 





text( -60, 10, sprintf('$U = %.0f$', U), 'interpreter', 'Latex', 'FontSize', fsize);

text( -70, 45, sprintf('(a)'), 'interpreter', 'Latex', 'FontSize', fsize);

c=colorbar('southoutside', 'Direction','reverse');
%colorbar('Direction','reverse')
c.Label.String = '$\delta n(x,t)/\mu$';
c.Label.Interpreter = 'latex';
c.Label.FontSize = fsize+5;


fname = sprintf('average_occupation_2D_U_%.2f_mu_%.4f.pdf',  U,  mu);

set(gcf,'paperunits','in');
set(gcf,'papersize',[7.3,5.3]); % Desired outer dimensions
% of figure

hfig = gcf;
print(hfig,'-bestfit','-dpdf',fname);




%% Generates the 2d density plot for the current for U=0, ballistic.

mu = 0.0005;
L = 200;
Mmax = 200;
U = 0.00;
dtm = 0.2;

datafile = ['./TEBD_no_leads_L_', num2str(L),'_U_',num2str(U, '%.2f'), '_mu_', num2str(mu, '%.4f'), '_Mmax_', num2str(Mmax) , '/datafile.mat']
load(datafile);

average_current_TEBD = (average_current_TEBD)/(mu);
h2=figure('color','white','units','inches','position',[1 1 7 5]);
hold on; 

size_x = size(average_current_TEBD, 2);
size_t = size(average_current_TEBD, 1);

sites = [-size_x/2:size_x/2-1];
times = [1:size_t]*dtm;
pcolor(sites, times, average_current_TEBD);

shading interp;

xlim ([-78,78]);
ylim([0,50]);

fsize = 20;
set(gca, 'Fontsize', fsize);
xlabel(sprintf('$x$'),'Interpreter', 'latex', 'FontSize',fsize+5);
ylabel(sprintf('$t$'),'Interpreter', 'latex', 'FontSize',fsize+5);


x = [0:-0.1:-70];
y = -x;
plot(x,y,'Color','red','LineStyle','--')

x = [0:0.1:70];
y = x;
line(x,y,'Color','red','LineStyle','--')
a=text( -60, 10, sprintf('$U = %.0f$', U), 'interpreter', 'Latex', 'FontSize', fsize, 'Color', 'w');

text( -70, 45, sprintf('(c)'), 'interpreter', 'Latex', 'FontSize', fsize, 'Color', 'w');

c=colorbar('southoutside', 'Direction','reverse');
%colorbar('Direction','reverse')
c.Label.String = '$j(x,t)/\mu$';
c.Label.Interpreter = 'latex';
c.Label.FontSize = fsize+5;


fname = sprintf('average_current_2D_U_%.2f_mu_%.4f.pdf',  U,  mu);

set(gcf,'paperunits','in');
set(gcf,'papersize',[7.3,5.3]); % Desired outer dimensions
% of figure

hfig = gcf;
print(hfig,'-bestfit','-dpdf',fname);



%% Generates the 2d density plot for density average for U\ne 0, superdiffusive
close all;
clear all;

mu = 0.0005;
L = 160;
Mmax = 200;
U = 1.00;
dtm = 0.2;

datafile = ['./TEBD_no_leads_L_', num2str(L),'_U_',num2str(U, '%.2f'), '_mu_', num2str(mu, '%.4f'), '_Mmax_', num2str(Mmax) , '/datafile.mat']
load(datafile);

average_occupation_TEBD = (average_occupation_TEBD-1)/(mu/2);
h2=figure('color','white','units','inches','position',[1 1 7 5]);
hold on; 

size_x = size(average_occupation_TEBD, 2);
size_t = size(average_occupation_TEBD, 1);

sites = [-size_x/2:size_x/2-1];
times = [1:size_t]*dtm;
pcolor(sites, times, average_occupation_TEBD);

shading interp;

xlim ([-78,78]);
ylim([0,50]);

fsize = 20;
set(gca, 'Fontsize', fsize);
xlabel(sprintf('$x$'),'Interpreter', 'latex', 'FontSize',fsize+5);
ylabel(sprintf('$t$'),'Interpreter', 'latex', 'FontSize',fsize+5);


x = [0:0.1:70];
y = x.^1.5;
plot(-x,y,'Color','k','LineStyle','--')

x = [0:0.1:70];
y = x.^1.5;
line(x,y,'Color','k','LineStyle','--')

n1_half = [];
n2_half = [];
for t=5:size_t
    [x, idx] = unique(average_occupation_TEBD(t, :));
    
    y = sites(idx);

    x_half = interp1(x,y, 0.5);
    n1_half =[n1_half; x_half, times(t)]; 

    x_half = interp1(x, y, -0.5);
    n2_half =[n2_half; x_half, times(t)]; 
end

plot(n1_half(1:10:end,1), smooth(n1_half(1:10:end,2)), 'ro' );
plot(n2_half(1:10:end,1)+1, n2_half(1:10:end,2), 'ro' );



c=colorbar('southoutside', 'Direction','reverse');
%colorbar('Direction','reverse')
c.Label.String = '$\delta n(x,t)/\mu$';
c.Label.Interpreter = 'latex';
c.Label.FontSize = fsize+5;

text( -60, 10, sprintf('$U = %.1f$', U), 'interpreter', 'Latex', 'FontSize', fsize);
text( -70, 45, sprintf('(b)'), 'interpreter', 'Latex', 'FontSize', fsize);

fname = sprintf('average_occupation_2D_U_%.2f_mu_%.4f.pdf',  U,  mu);

set(gcf,'paperunits','in');
set(gcf,'papersize',[7.3,5.3]); % Desired outer dimensions
% of figure

hfig = gcf;
print(hfig,'-bestfit','-dpdf',fname);





%% Current across the interface as function of time for U=0; 

close all;
clear all;
mu = 0.0005;
L = 200;
Mmax = 200;
U = 0.00;
dtm = 0.2;
step =1;
datafile = ['./TEBD_no_leads_L_', num2str(L),'_U_',num2str(U, '%.2f'), '_mu_', num2str(mu, '%.4f'), '_Mmax_', num2str(Mmax) , '/datafile.mat']
load(datafile);

size_t = size(average_current_TEBD, 1);
times = [1:size_t]'*dtm;

current = average_current_TEBD(:,L/2)/mu;

h2=figure('color','white','units','inches','position',[1 1 7 5]);
hold on;

plot(times(1:step:end), current(1:step:end),...
    'linestyle','None',...
    'Marker', 'o',...
    'MarkerSize',6,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor','w', 'DisplayName', 'Numerics');
fsize = 20;


set(gca, 'Fontsize', fsize);

current = average_current_TEBD(:,L/2)/mu;
plot(times(1:step:end), current(1:step:end),'b-', 'DisplayName', 'Analytical');

xlim([0, 60]);

ylabel(sprintf('$\\langle  j(x=0,t)/\\mu \\rangle $'),'Interpreter', 'latex', 'FontSize',fsize+3);
xlabel(sprintf('$t$'),'Interpreter', 'latex', 'FontSize',fsize+3)

leg=legend;
leg.Location = 'SouthEast';
leg.Position= [0.55 0.7 0.25 0.15]
box on;

text( 10, 0.8, sprintf('$U = %.1f$', U), 'interpreter', 'Latex', 'FontSize', fsize);

text( 5, 0.1, sprintf('(a)', U), 'interpreter', 'Latex', 'FontSize', fsize);

axes('Position',[.6 .32 .35 .3])
hold on; 
box on

step = 2;

cumulative_current = cumtrapz(times, current);
plot(times(1:step:end), cumulative_current(1:step:end),...
    'linestyle','None',...
    'Marker', 'o',...
    'MarkerSize',6,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor','w', 'DisplayName', 'Numerics');

fig = gca;


ylabel(sprintf('$\\langle \\Delta j(0,t)/\\mu \\rangle $'),'Interpreter', 'latex', 'FontSize',fsize+3);
xlabel(sprintf('$t$'),'Interpreter', 'latex', 'FontSize',fsize+3)


% fitting the exponent. 
t_start = 200; 
t_end = 500; 
x = times(t_start:t_end);
y = cumulative_current(t_start:t_end);
pol = polyfit(log(x),log(y), 1);
m = pol(1);
k = pol(2);

cumulative_current_fit = times.^m.*exp(k);

plot(times, cumulative_current_fit, 'k-', 'DisplayName', 'Fit');

fig.XScale ='log';
fig.YScale = 'log';
fig.FontSize = fsize;
fig.XTick=[1, 10, 100];
fig.YTick = [1];



text( 8, 0.8, sprintf('$\\alpha = %.1f$', m), 'interpreter', 'Latex', 'FontSize', fsize);



fname = sprintf('current_U_%.2f_mu_%.4f.pdf',  U,  mu);
set(gcf,'paperunits','in');
set(gcf,'papersize',[7.3,5.3]); % Desired outer dimensions
% of figure
leg1 =legend;
leg1.Position=[0.58 0.48 0.1 0.1];
hfig = gcf;
print(hfig,'-bestfit','-dpdf',fname);










%% Current across the interface as function of time for U\ne 0; 

close all;
mu = 0.0005;
L = 400;
Mmax = 200;
U = 1.00;
dtm = 0.2;
step =1;
datafile = ['./TEBD_no_leads_L_', num2str(L),'_U_',num2str(U, '%.2f'), '_mu_', num2str(mu, '%.4f'), '_Mmax_', num2str(Mmax) , '/datafile.mat']
load(datafile);

size_t = size(average_current_TEBD, 1);
times = [1:size_t]*dtm;
current = average_current_TEBD/mu;

h2=figure('color','white','units','inches','position',[1 1 7 5]);
hold on;

plot(times(1:step:end), current(1:step:end),...
    'linestyle','None',...
    'Marker', 'o',...
    'MarkerSize',6,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor','w', 'DisplayName', 'Numerics');
fsize = 20;

fig = gca;
fig.XScale = 'Log';
fig.YScale = 'Log';


fig.FontSize = fsize;

current = average_current_TEBD/mu;

xlim([0.2, 220]);
ylim([0,1.5]);

ylabel(sprintf('$\\langle  j(0,t)/\\mu \\rangle $'),'Interpreter', 'latex', 'FontSize',fsize+3);
xlabel(sprintf('$t$'),'Interpreter', 'latex', 'FontSize',fsize)





% fitting the exponent. 
t_start = 100; 
t_end = 1000; 
x = times(t_start:t_end);
y = current(t_start:t_end);
pol = polyfit(log(x),log(y), 1);
m = pol(1)
k = pol(2);

time_interval = 10:0.1:400
current_fit = time_interval.^m.*exp(k);

plot(time_interval, current_fit, 'k-', 'DisplayName', 'Fit $\langle j(0,t)\rangle \sim t^{-1/3}$');

leg=legend;
leg.Interpreter='latex';
 leg.Location = 'SouthEast';
 leg.Position= [0.4 0.2 0.25 0.15]


 box on;

text( 0.25, 1, sprintf('$U = %.1f$', U), 'interpreter', 'Latex', 'FontSize', fsize+3);
text( 0.3, 0.14, sprintf('(b)'), 'interpreter', 'Latex', 'FontSize', fsize);


axes('Position',[.59 .6 .3 .3])
hold on; 
box on

step = 10;

cumulative_current = cumtrapz(times, current);
plot(times(1:step:end), cumulative_current(1:step:end),...
    'linestyle','None',...
    'Marker', 'o',...
    'MarkerSize',6,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor','w', 'DisplayName', 'Numerics');

fig = gca;


ylabel(sprintf('$\\langle \\Delta j(0,t)/\\mu \\rangle $'),'Interpreter', 'latex', 'FontSize',fsize+3);
xlabel(sprintf('$t$'),'Interpreter', 'latex', 'FontSize',fsize+3)


% fitting the exponent. 
t_start = 700; 
t_end = 1000; 
x = times(t_start:t_end);
y = cumulative_current(t_start:t_end);
pol = polyfit(log(x),log(y), 1);
m = pol(1);
k = pol(2);

cumulative_current_fit = times.^m.*exp(k);

plot(times, cumulative_current_fit, 'k-', 'DisplayName', 'Fit');

fig.XScale ='log';
fig.YScale = 'log';
fig.FontSize = fsize;
fig.XTick=[1, 10, 100];
fig.YTick = [1];



text( 8, 0.7, sprintf('$\\alpha = %.2f$', m), 'interpreter', 'Latex', 'FontSize', fsize);



fname = sprintf('current_U_%.2f_mu_%.4f.pdf',  U,  mu);
set(gcf,'paperunits','in');
set(gcf,'papersize',[7.3,5.3]); % Desired outer dimensions
% of figure
leg1 =legend;
leg1.Position=[0.62 0.8 0.1 0.1];
hfig = gcf;
print(hfig,'-bestfit','-dpdf',fname);








%% Scaling exponent for U=0 and finite U; 

% U=0 case;
close all;
mu = 0.0005;
L = 160;
Mmax = 200;
U = 0.00;
dtm = 0.2;
step =1;
datafile = ['./TEBD_no_leads_L_', num2str(L),'_U_',num2str(U, '%.2f'), '_mu_', num2str(mu, '%.4f'), '_Mmax_', num2str(Mmax) , '/datafile.mat']
load(datafile);

size_t = size(average_current_TEBD, 1);
times = [1:size_t]'*dtm;
current = average_current_TEBD/mu;
cumulative_current = cumtrapz(times, current);

 exponents0 = diff(log(cumulative_current))./diff(log(times));


h2=figure('color','white','units','inches','position',[1 1 7 5]);
hold on;

p1=plot(times(1:step:end-1), exponents0(1:step:end),...
    'linestyle','None',...
    'Marker', 'o',...
    'MarkerSize',6,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor','w');
fsize = 16;



fig = gca; 
fig.FontSize = fsize;

current = average_current_TEBD/mu;
p2=plot(times(1:step:end-1), exponents0(1:step:end),'k-');



ylabel(sprintf('$\\alpha(t)  $'),'Interpreter', 'latex', 'FontSize',fsize+3);
xlabel(sprintf('$t$'),'Interpreter', 'latex', 'FontSize',fsize+3)


 
 box on;

text( 10, 2.5, sprintf('$U = %.1f$', U), 'interpreter', 'Latex', 'FontSize', fsize);



% U = 1.0 case;



mu = 0.0005;
L = 400;
Mmax = 200;
U = 1.00;
dtm = 0.2;
step1 =5;
datafile = ['./TEBD_no_leads_L_', num2str(L),'_U_',num2str(U, '%.2f'), '_mu_', num2str(mu, '%.4f'), '_Mmax_', num2str(Mmax) , '/datafile.mat']
load(datafile);



size_t = size(average_current_TEBD, 1);
times = [1:size_t]'*dtm;
current = average_current_TEBD/mu;
cumulative_current = cumtrapz(times, current);

 exponents0 = diff(log(cumulative_current))./diff(log(times));


p3=plot(times(1:step1:end-500), exponents0(1:step1:end-499),...
    'linestyle','None',...
    'Marker', 'd',...
    'MarkerSize',6,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor','w');

ylim([0.5:1.5]);
xlim([0, 300]);


%draw some lines to guide the eye. 

x = [0 300];
y = [1 1];
line(x,y,'Color','red','LineStyle','--');

x = [0 300];
y = [0.66 0.66];
line(x,y,'Color','blue','LineStyle','--');


leg = legend([p1, p2, p3], 'U = 0, numerics','U = 0, analytical','U = 1.0, numerics',''); 


text( 10, 2.5, sprintf('$U = %.1f$', U), 'interpreter', 'Latex', 'FontSize', fsize);

text(170, 1.03, sprintf('$\\alpha_b = %.1f$', 1.0), 'interpreter', 'Latex', 'FontSize', fsize);
text(240, 0.7, sprintf('$\\alpha_{sd} = {2\\over 3}$'), 'interpreter', 'Latex', 'FontSize', fsize);

% text(50, 1.05, sprintf('$U = %.1f$', 0), 'interpreter', 'Latex', 'FontSize', fsize);
% text(80, 0.75, sprintf('$U = %.1f$', 1.0), 'interpreter', 'Latex', 'FontSize', fsize);

fname = sprintf('scaling_exponent_mu_%.4f.pdf',   mu);
set(gcf,'paperunits','in');
set(gcf,'papersize',[7.3,5.3]); % Desired outer dimensions
% of figure
hfig = gcf;
print(hfig,'-bestfit','-dpdf',fname);









%% Generates the scaling for  density average for U=0, ballistic.
close all;
mu = 0.0005;
L = 200;
Mmax = 200;
U = 0.00;
dtm = 0.2;

datafile = ['./TEBD_no_leads_L_', num2str(L),'_U_',num2str(U, '%.2f'), '_mu_', num2str(mu, '%.4f'), '_Mmax_', num2str(Mmax) , '/datafile.mat']
load(datafile);

average_occupation_TEBD = (average_occupation_TEBD-1)/(mu/2);
h2=figure('color','white','units','inches','position',[1 1 7 5]);
hold on; 

size_x = size(average_occupation_TEBD, 2);
size_t = size(average_occupation_TEBD, 1);

sites = [-size_x/2:size_x/2-1];
times = [1:size_t]*dtm;

times_idx = [ 100,200 300 ];

symbols=['s','d','o','v','s','d','o','v', 's','d','o','v'];
colors = ['k', 'b','r','m','r', 'b','k','m','r', 'b','k','m'];

for j=1:numel(times_idx)
 
    plot(sites/times(times_idx(j)), average_occupation_TEBD (times_idx(j), :),...
        'linestyle','None',...
        'Marker', symbols(j),...
        'MarkerSize', 5,...
        'MarkerEdgeColor',colors(j),...
        'MarkerFaceColor','none', 'DisplayName', sprintf('t = %.1f', times(times_idx(j))));


end

xlim ([-3,3]);
ylim([-1,1]);

fsize = 20;
set(gca, 'Fontsize', fsize);
xlabel(sprintf('$x/t$'),'Interpreter', 'latex', 'FontSize',fsize+5);
ylabel(sprintf('$ \\delta n(x,t)/\\mu$'),'Interpreter', 'latex', 'FontSize',fsize+5);


leg=legend;
leg.Position= [0.6 0.7 0.2 0.15]
box on;

text(-2, -0.5, sprintf('$U = %.1f$', U), 'interpreter', 'Latex', 'FontSize', fsize);
text(-2.7, 0.8, sprintf('(a)'), 'interpreter', 'Latex', 'FontSize', fsize);


fname = sprintf('scaling_occupation_U_%.2f_mu_%.4f.pdf',  U,  mu);

set(gcf,'paperunits','in');
set(gcf,'papersize',[7.3,5.3]); % Desired outer dimensions
% of figure

hfig = gcf;
print(hfig,'-bestfit','-dpdf',fname);



%% Generates the scaling for  current average for U=0, ballistic.
close all;
mu = 0.0005;
L = 200;
Mmax = 200;
U = 0.00;
dtm = 0.2;

datafile = ['./TEBD_no_leads_L_', num2str(L),'_U_',num2str(U, '%.2f'), '_mu_', num2str(mu, '%.4f'), '_Mmax_', num2str(Mmax) , '/datafile.mat']
load(datafile);

average_current_TEBD = (average_current_TEBD)/mu;
h2=figure('color','white','units','inches','position',[1 1 7 5]);
hold on; 

size_x = size(average_occupation_TEBD, 2);
size_t = size(average_occupation_TEBD, 1);

sites = [-size_x/2+1:size_x/2-1];
times = [1:size_t]*dtm;

times_idx = [ 100,200 300 ];

symbols=['s','d','o','v','s','d','o','v', 's','d','o','v'];
colors = ['k', 'b','r','m','r', 'b','k','m','r', 'b','k','m'];

for j=1:numel(times_idx)
 
    plot(sites/times(times_idx(j)), average_current_TEBD (times_idx(j), :),...
        'linestyle','None',...
        'Marker', symbols(j),...
        'MarkerSize', 5,...
        'MarkerEdgeColor',colors(j),...
        'MarkerFaceColor','none', 'DisplayName', sprintf('t = %.1f', times(times_idx(j))));


end

xlim ([-3,3]);
ylim([0,0.9]);

fsize = 20;
set(gca, 'Fontsize', fsize);
xlabel(sprintf('$x/t$'),'Interpreter', 'latex', 'FontSize',fsize+5);
ylabel(sprintf('$ \\langle j(x,t)\\rangle /\\mu$'),'Interpreter', 'latex', 'FontSize',fsize+5);


leg=legend;
leg.Position= [0.65 0.7 0.2 0.15]
box on;

text(-2.5, 0.4, sprintf('$U = %.1f$', U), 'interpreter', 'Latex', 'FontSize', fsize);
text(-2.7, 0.8, sprintf('(c)'), 'interpreter', 'Latex', 'FontSize', fsize);


fname = sprintf('scaling_current_U_%.2f_mu_%.4f.pdf',  U,  mu);

set(gcf,'paperunits','in');
set(gcf,'papersize',[7.3,5.3]); % Desired outer dimensions
% of figure

hfig = gcf;
print(hfig,'-bestfit','-dpdf',fname);


%% Generates the scaling for  density average for U\ne 0, superdiffusive.
close all;
mu = 0.0005;
L = 400;
Mmax = 200;
U = 1.00;
dtm = 0.2;

datafile = ['./TEBD_no_leads_L_', num2str(L),'_U_',num2str(U, '%.2f'), '_mu_', num2str(mu, '%.4f'), '_Mmax_', num2str(Mmax) , '/datafile.mat']
load(datafile);

average_occupation_TEBD = (average_occupation_TEBD-1)/(mu/2);
h2=figure('color','white','units','inches','position',[1 1 7 5]);
hold on; 

size_x = size(average_occupation_TEBD, 2);
size_t = size(average_occupation_TEBD, 1);

sites = [-size_x/2:size_x/2-1];
times = [1:size_t]*dtm;

times_idx = [ 400,800 1000 1200];

symbols=['s','d','o','v','s','d','o','v', 's','d','o','v'];
colors = ['k', 'b','r','g','r', 'b','k','m','r', 'b','k','m'];

for j=1:numel(times_idx)
 
    plot(sites/times(times_idx(j)).^0.66, average_occupation_TEBD (times_idx(j), :),...
        'linestyle','None',...
        'Marker', symbols(j),...
        'MarkerSize', 5,...
        'MarkerEdgeColor',colors(j),...
        'MarkerFaceColor','none', 'DisplayName', sprintf('t = %.1f', times(times_idx(j))));


end

xlim ([-8,8]);
ylim([-1,1]);

fsize = 20;
set(gca, 'Fontsize', fsize);
xlabel(sprintf('$x/t^{2/3}$'),'Interpreter', 'latex', 'FontSize',fsize+5);
ylabel(sprintf('$ \\delta n(x,t)/\\mu$'),'Interpreter', 'latex', 'FontSize',fsize+5);


leg=legend;
leg.Position= [0.6 0.7 0.2 0.15]
box on;

text(-5, -0.5, sprintf('$U = %.1f$', U), 'interpreter', 'Latex', 'FontSize', fsize);
text(-7.5, 0.8, sprintf('(b)', U), 'interpreter', 'Latex', 'FontSize', fsize);


fname = sprintf('scaling_occupation_U_%.2f_mu_%.4f.pdf',  U,  mu);

set(gcf,'paperunits','in');
set(gcf,'papersize',[7.3,5.3]); % Desired outer dimensions
% of figure

hfig = gcf;
print(hfig,'-bestfit','-dpdf',fname);





%% Gradient of the density to check for the scaling. 


 


close all;
clear all;

a= 0.67;
b = 0.67;

mu = 0.0005;
L = 400;
Mmax = 200;
U = 1.00;
dtm = 0.2;

datafile = ['./TEBD_no_leads_L_', num2str(L),'_U_',num2str(U, '%.2f'), '_mu_', num2str(mu, '%.4f'), '_Mmax_', num2str(Mmax) , '/datafile.mat']
load(datafile);

average_occupation_TEBD = (average_occupation_TEBD-1)/(mu/2);

h2=figure('color','white','units','inches','position',[1 1 7 5]);
hold on; 

size_x = size(average_occupation_TEBD, 2);
size_t = size(average_occupation_TEBD, 1);

sites = [-size_x/2+1:size_x/2-1];
times = [1:size_t]*dtm;

times_idx = [ 1000 1200 1400 ];

symbols=['s','d','o','v','s','d','o','v', 's','d','o','v'];
colors = ['m', 'b','r','g','r', 'b','k','m','r', 'b','k','m'];

for j=1:numel(times_idx)
 
    plot(b*sites(1:2:end-1)/times(times_idx(j)).^(2/3), -times(times_idx(j)).^(2/3)*diff(average_occupation_TEBD (times_idx(j), 1:2:end)),...
        'linestyle','None',...
        'Marker', symbols(j),...
        'MarkerSize', 5,...
        'MarkerEdgeColor',colors(j),...
        'MarkerFaceColor','none', 'DisplayName', sprintf('t = %.1f', times(times_idx(j))));


end

%xlim ([-8,8]);
%ylim([-1,1]);

fig =gca; 
%fig.YScale='log';

% plot a gaussian 
x= -5:0.01:5;
gauss = @(x,mu,sig,amp,vo) amp*exp(-(((x-mu).^2)/(2*sig.^2)))+vo;

y = gauss(x, 0.0, 1, 1.2 ,0 );
plot(x,y, 'k-','Linewidth',1.3, 'DisplayName', 'Gaussian fit');




fsize = 16;
set(gca, 'Fontsize', fsize);
xlabel(sprintf('$b x/t^{2/3}$'),'Interpreter', 'latex', 'FontSize',fsize+5);
ylabel(sprintf('$ t^{2/3}\\Delta n(x,t)/\\mu$'),'Interpreter', 'latex', 'FontSize',fsize+5);


leg=legend;
leg.Position= [0.2 0.7 0.2 0.15];
leg.FontSize = 16;      
box on;


fig.XLim=[-4,7];
fig.YLim=[0, 2];
text(-3.5, 0.75, sprintf('$U = %.1f$', U), 'interpreter', 'Latex', 'FontSize', fsize);
%text(-7.5, 0.8, sprintf('(b)', U), 'interpreter', 'Latex', 'FontSize', fsize);



axes('Position',[.58 .4 .3 .5])
hold on; 
box on

fig1 = gca;

for j=1:numel(times_idx)
 
    plot(b*sites(1:2:end-1)/times(times_idx(j)).^(2/3), -times(times_idx(j)).^(2/3)*diff(average_occupation_TEBD (times_idx(j), 1:2:end)),...
        'linestyle','None',...
        'Marker', symbols(j),...
        'MarkerSize', 5,...
        'MarkerEdgeColor',colors(j),...
        'MarkerFaceColor','none', 'DisplayName', sprintf('t = %.1f', times(times_idx(j))));


end
plot(x,y, 'k-','Linewidth',1.3, 'DisplayName', 'Gaussian fit');



xlabel(sprintf('$b x/t^{2/3}$'),'Interpreter', 'latex', 'FontSize',fsize);
ylabel(sprintf('$ t^{2/3}\\Delta n(x,t)/\\mu$'),'Interpreter', 'latex', 'FontSize',fsize);


fig1.YScale ='log'; 
fig1.XLim=[-3,3];

fname = sprintf('scaling_function_U_%.2f_mu_%.4f.pdf',  U,  mu);

set(gcf,'paperunits','in');
set(gcf,'papersize',[7.3,5.3]); % Desired outer dimensions
% of figure

hfig = gcf;
print(hfig,'-bestfit','-dpdf',fname);























