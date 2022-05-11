clc
clear all

N = 1000;
hami = sparse(N,N);

bound = 6;
x = linspace(-bound, bound,N);
dx = x(2) - x(1);
kin = sparse(N,N);

m = 1;
hbar = 1;

mtx = 2 * eye(N);
for i = 2:N
    mtx(i-1,i) = -1;
    mtx(i,i-1) = -1;
end

kin = hbar^2 / (2 * m * dx^2) * mtx;

pot = sparse(N,N);
a = -1;
U =  0.5 * a * x.^2 + 0.25 * x.^4;
U = U - min(U);

figure(1)
clf(figure(1))
hold on
plot(x,U)
hold off

for i = 1:N
   pot(i,i) = U(i);
end
%pot = pot - min(pot);

hami = kin + pot;
[psi,E] = eig(hami);

EE = diag(E);

figure(3)
hold on
plot(EE)
hold off

figure(2)
clf(figure(2))
hold on
%legend(['a','b'],'Location','southwest')

plot(x,psi(:,1))
plot(x,psi(:,2))

hold off
legend(num2str(EE(1)),num2str(EE(2)))



k = 2 *( EE(1) - U);
k2 = 2* (EE(2) -U);

P = zeros(N,1);
P(1) = 0;
P(2) = 0.0000001;

P2 = zeros(N,1);
P2(1) = P(1);
P2(2) = P(2);

for i = 3:N
    P(i) = (-k(i)*dx*dx + 2) * P(i-1) - P(i-2);
    P2(i) = (-k2(i)*dx*dx + 2) * P2(i-1) - P2(i-2);
end
    %%
    
figure(4)
clf(figure(4))
hold on
plot(x,P)
plot(x,P2 + EE(2))
hold off
   
    