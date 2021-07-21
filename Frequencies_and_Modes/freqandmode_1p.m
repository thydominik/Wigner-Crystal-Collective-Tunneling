clc
clear all

%constants
l_d = 161.07;                  %length scale of the potential
E_p = 0.471738;% * 1.602176487*10^-19 * 10^-3;                %potential scaling factor
alpha = 6.343306112300098;                  %x^2 parameter
khi_0 = sqrt(alpha);            %equlibrium position
me = 9.1093837015*10^-31;
m = 0.0062 * me;
hbar = 1.0545718  * 10^-34;

Vdd = E_p/2 * (3*khi_0 + alpha);  %second derivative of the potential.

w = sqrt(Vdd);
disp(['\omega = ',num2str(w),' from experiment data'])


N = 1000;
a = linspace(-20,20,N);

figure(2)
hold on
clf(figure(2))
plot(a,sqrt(-(a)))

hold off




for i=1:N
    w(i) = sqrt((1 * (3*sqrt(abs(a(i)))^2 + a(i))));
end




%mivel E_p nek van egyedül a potenciálban jelenleg mértékegysége emiatt
%omega meV ban van megadva!

% gamow factor calculatio

x = linspace(-khi_0*l_d,khi_0*l_d,100000);
V = E_p * ((1/4) * (x/l_d).^4 - (alpha/2) * (x/l_d).^2);

E_0 = hbar * w(500)/(1.602176487*10^-19);

g = sqrt(2*m) .* sqrt((E_0 - V));
% figure(2)
% clf(figure(2))
% hold on
% plot(x,g)
% hold off

G = 0;
for i = 1:length(x)
    G = G + g(i);
end
format long
e_g = exp(G);
disp(['G = ',num2str(G),';  exp Gamow factor = ',num2str(e_g)])


figure(1)
clf(figure(1))
hold on
plot(-a,w,'r','LIneWidth',2)
xlabel('$$\tilde{a}$$', 'Interpreter', 'LaTeX', 'FontSize', 16)
ylabel('$\hbar\omega [meV]$', 'Interpreter', 'latex','FontSize',16)
%title('\omega^2 (\alpha)','FontSize',17)
hold off
    