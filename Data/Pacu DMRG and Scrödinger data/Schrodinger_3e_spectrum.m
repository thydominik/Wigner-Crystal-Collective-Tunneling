close all;
clear all;
Ne=3;
eta = 20.0;
beta = 0.3;
levels = 11;
degeneracy = 8; 

xMin = -6; % xmin in units of ld
xMax = 6; %xmax in units of ld;
N=100;
Norb = 8; 

x=linspace(xMin, xMax, N);

filename = sprintf('E_Schrodinger_%de_eta_%.2f_N_%d_beta_%.3f.dat', Ne, eta, N, beta);

E= load(filename);
alpha =E(:,1);

symbol1 = ['+','+', '+', '+'];
color1 = ['r','r', 'r', 'r'];

max_color = 4;


h1=figure('color','white','units','inches','position',[1 1 7 4]);hold on;
r= 12;

for j=2:degeneracy:levels
    
    plot(-alpha, E(:,j),...
        'linestyle','-',...
        'LineWidth', 1);
    
end

fsize = 20;
%set(gca, 'YScale', 'log')
xlabel('$ -\alpha $','Interpreter', 'latex', 'FontSize',fsize+3);
ylabel(sprintf('$E_n (\\eta = %.1f) $', eta),'Interpreter', 'latex', 'FontSize',fsize+3);


if 1
    set(gcf, 'PaperPositionMode', 'auto');
    fname = sprintf('Energy_spectrum_eta_%.2f_beta_%.2f.pdf', eta, beta);
    
    set(gcf,'paperunits','in');
    set(gcf,'papersize',[7.3,4.3]); % Desired outer dimensions
    % of figure
    
    hfig = gcf;
    print(hfig,'-dpdf',fname);
    
end



h2=figure('color','white','units','inches','position',[1 1 7 4]);hold on;
r= 12;
hold on; 

load('Delta_E_DMRG_Norb_8_eta_20.00.mat');
    
plot(-alpha_list ,delta_E_DMRG,  'o','LineWidth',1,...
        'MarkerEdgeColor','r',...
        'MarkerFaceColor','w',...
        'MarkerSize',7);   
    

 
load ('Splitting_instanton_eta_20.mat');
plot(-delta(:,1) ,delta(:,2),  's-','LineWidth',1,...
        'MarkerEdgeColor','b',...
        'MarkerFaceColor','w',...
        'MarkerSize',7);  

plot(-alpha, E(:,degeneracy+2)-E(:,2),...
        'Color', 'k', 'linestyle','-',...
        'LineWidth', 2);
        
fsize = 20;
ylim([-1, 8])
set(gca, 'YScale', 'log')
%set(gca, 'xScale', 'log')

xlabel('$ -\alpha $','Interpreter', 'latex', 'FontSize',fsize+3);
ylabel(sprintf('$\\Delta E (\\eta = %.1f) $', eta),'Interpreter', 'latex', 'FontSize',fsize+3);


   leg{1}=sprintf('$DMRG, N_{orb}=%d$', Norb);
    leg{2}=sprintf('$Instanton$');
    leg{3} = sprintf('$ED\\, \\, Schrodinger$')
    hleg=legend(leg,  'Location', 'SouthWest');
    legend('boxoff');
    set(hleg,'Interpreter','latex', 'FontSize', fsize);
   

if 1
    set(gcf, 'PaperPositionMode', 'auto');
    fname = sprintf('Energy_splitting_eta_%.2f_beta_%.2f.pdf', eta, beta);
    
    set(gcf,'paperunits','in');
    set(gcf,'papersize',[7.3,4.3]); % Desired outer dimensions
    % of figure
    
    hfig = gcf;
    print(hfig,'-dpdf',fname);
    
end


return; 

alpha_values = size(E, 1);



for k = 1:5:alpha_values
    h2=figure('color','white','units','inches','position',[1 1 7 4]);hold on;
    
    alpha_U = E(k, 1);
    
    U= 0.25*x.^4 + 0.5*alpha_U*x.^2;
    U = U-min(U);
    
    plot(x, U);

    x1 = linspace(xMin/2, xMax/2, N);
    
    for j=2:degeneracy:levels
        plot(x1, E(k, j)*ones(1, N));
    end
    
    xlim([-6,6]);
    ylim([0, 50]);
    %set(gca, 'YScale', 'log');
    xlabel('$ x$','Interpreter', 'latex', 'FontSize',fsize+3);
    ylabel(sprintf('$ E_n (\\eta = %.1f, \\alpha=%.1f) $', eta, alpha_U),'Interpreter', 'latex', 'FontSize',fsize+3);

    if 1
    set(gcf, 'PaperPositionMode', 'auto');
    fname = sprintf('Energy_spectrum_alpha_%.1f_eta_%.2f_beta_%.2f.pdf',alpha_U, eta, beta);
    
    set(gcf,'paperunits','in');
    set(gcf,'papersize',[7.3,4.3]); % Desired outer dimensions
    % of figure
    
    hfig = gcf;
    print(hfig,'-dpdf',fname);
    
   end


    
    
end
    


