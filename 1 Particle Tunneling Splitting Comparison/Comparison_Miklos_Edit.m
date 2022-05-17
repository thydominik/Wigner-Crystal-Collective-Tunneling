clc
clear all

%1 particle tunneling calculations for a feq methods
%Methods:
% Exact with Schrödinger eq.
% Landau's formula
% Colemans formula
% Instanton Claulcation
% Milnikov's formula

N       = 200;
Nalpha  = 100;
alpha   = linspace(0, 10, Nalpha);
Action = 2*sqrt(2)/3 * alpha.^(1.5);


m       = 1;
hbar    = 1;
 
splits = zeros(Nalpha, 5);
%First the Landau method:--------------------------------------------------
for i = 1:Nalpha
    

    x       = linspace(-sqrt(alpha(i)), sqrt(alpha(i)), N);
    V       = 0.25 * (alpha(i) - x.^2).^2;
    p       = abs(sqrt(2 * m * V));
    omega   =  sqrt(2 * abs(alpha(i))); %sqrt((1 * (3*sqrt(abs(alpha(i)))^2 - alpha(i))));
    
    int_p = 0;
    for j = 1:length(x)- 1
        int_p = int_p + (x(j+1) - x(j)) * (p(j) + p(j+1))/2;
    end

    splits(i, 1) = omega * hbar / pi * exp(-1/hbar * int_p);
end

%exact method:-------------------------------------------------------------
bound = 6;
x = linspace(-bound, bound,N);
dx = x(2) - x(1);
kin = sparse(N,N);

mtx = 2 * eye(N);
for i = 2:N
    mtx(i-1,i) = -1;
    mtx(i,i-1) = -1;
end

kin = hbar^2 / (2 * m * dx^2) * mtx;

for i = 1:Nalpha
   U = 0.25 * (alpha(i) - x.^2).^2;
   pot = diag(U);
   hami = kin + pot;
   [psi,E] = eig(hami);
   EE(i,:) = diag(E);
end

for i = 1:Nalpha
   splits(i, 2) = EE(i,2) - EE(i,1);    
end

%Coleman:------------------------------------------------------------------
% Action = 0.924 * alpha.^(1.5);
% for i = 1:Nalpha
%     xint = linspace(0,alpha(i),N);
%     dxp2 = (xint(5)-xint(4))/2;
%     
%     omega   =  sqrt(2 * abs(alpha(i))); %sqrt((1 * (3*sqrt(abs(alpha(i)))^2 - alpha(i))));
%     V = 0.25 * ( alpha(i) - xint.^2).^2;
%     A = zeros(Nalpha,1);
%     
%     func1 = (m*omega)./sqrt(2*m*V);
%     func2 = -1./(alpha(i) - xint);
%         
%     for j = 2:N-2
%         A(i) = A(i) + (func1(j)+func1(j+1)) * dxp2;
%         A(i) = A(i) + (func2(j)+func2(j+1)) * dxp2;
%     end  
%     
% end
% for k = 1:Nalpha
%    omega   =  sqrt(2 * abs(alpha(i))); %sqrt((1 * (3*sqrt(abs(alpha(i)))^2 - alpha(i))));
%    splits(k, 3) = hbar * omega * sqrt((m*omega*alpha(k)^2)/(pi*hbar)) * exp(A(k) - Action(k)/hbar);    
% end

% Instanton:---------------------------------------------------------------

for i = 1:Nalpha
    omega   =  sqrt(2 * abs(alpha(i))); %sqrt((1 * (3*sqrt(abs(alpha(i)))^2 - alpha(i))));
    %splits(i, 4) = 2 * hbar * sqrt(12) * omega * sqrt(Action(i)/(2 * pi * hbar)) * exp(- Action(i) / hbar);
    
    splits(i, 4) = 2 * hbar * sqrt(12) * omega * sqrt(Action(i)/(2 * pi * hbar)) * exp(- Action(i) / hbar);
end

% Milnikov:-----------------------------------------------------------------
y = linspace(-100, 0, N);
dy = y(2) - y(1);
chi = tanh(y /sqrt(2));
for k = 1:Nalpha
    V       = 0.25 * (x.^2 - alpha(k)).^2;
    p       = sqrt(2 * m * V);
    p0      = sqrt(2* m * (0.25 * (alpha(k))^2));
    omega   =  sqrt(2 * abs(alpha(k))); %sqrt((1 * (3*sqrt(abs(alpha(i)))^2 - alpha(i))));
    %splits(k, 5) = sqrt(omega/pi) * 4 * omega * sqrt(alpha(k)) * exp(-Action(k)/hbar); %* exp(A(k));
    splits(k, 5) = sqrt(omega/pi) * 4 * omega * sqrt(alpha(k)) * exp(-Action(k) /hbar); %* exp(A(k));
end
%%
figure(1)
clf(figure(1))
hold on
x = linspace(0,4,100);
func = 1.08 - 0.37 * x;
plot(x, func)
%plot(alpha, splits(:, 1),'.-', 'DisplayName', 'Landau')
plot(alpha, splits(:, 2),'.-', 'LineWidth', 2, 'DisplayName', 'Schrödinger equation')
%plot(alpha, splits(:, 3),'.-', 'DisplayName', 'Colemans approximation')
%plot(alpha, splits(:, 4),'.-', 'DisplayName', 'Instanton approximation')
plot(alpha, splits(:, 5),'.-', 'DisplayName', 'Milnikovs approximation')
legend
xlabel('| \alpha |')
ylabel('\Delta')
yline([10^-1 10^-2])
%set(gca, 'Xscale', 'log')
%set(gca, 'Yscale', 'log')
xlim([0 3])
ylim([0 1.1])
grid on

% figure(2)
% clf(figure(2))
% hold on
% plot(alpha, abs(splits(:, 2) - splits(:, 1)),'.-', 'DisplayName', 'Landau')
% plot(alpha, abs(splits(:, 2) - splits(:, 3)),'.-', 'DisplayName', 'Coleman')
% plot(alpha, abs(splits(:, 2) - splits(:, 4)),'.-', 'DisplayName', 'Instanton')
% plot(alpha, abs(splits(:, 2) - splits(:, 5)),'.-', 'DisplayName', 'Milnikov')
% legend
% set(gca, 'Yscale', 'log')
% hold off
% 
% figure(3)
% clf(figure(3))
% hold on
% plot(alpha(1:end-1), diff(splits(:, 1)),'.-', 'DisplayName', 'Landau')
% plot(alpha(1:end-1), diff(splits(:, 2)),'.-', 'LineWidth', 2, 'DisplayName', 'Schrödinger equation')
% plot(alpha(1:end-1), diff(splits(:, 3)),'.-', 'DisplayName', 'Colemans approximation')
% plot(alpha(1:end-1), diff(splits(:, 4)),'.-', 'DisplayName', 'Instanton approximation')
% plot(alpha(1:end-1), diff(splits(:, 5)),'.-', 'DisplayName', 'Milnikovs approximation')
% legend
% 
% hold off

%%
d = load('E_Schrodinger_3e_eta_20.00_beta_0.01_N_100.dat');
figure(1)
clf(figure(1))
hold on
D = abs(d(:, 2) - d(:, 3));
a = - d(:, 1);
plot(a, D)

hold off






