clc
clear all

N = 100;
hami = sparse(N,N);

bound = 5;
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
a = -3;
U =  0.5 * a * x.^2 + 0.25 * x.^4;
U = U - min(U);



for i = 1:N
   pot(i,i) = U(i);
end
%pot = pot - min(pot);

hami = kin + pot;
[psi,E] = eig(hami);

EE = diag(E);
splitting = EE(2) - EE(1);
filename = append('splitting_', num2str(abs(a)),'.txt');
save(filename, 'splitting','-ascii')