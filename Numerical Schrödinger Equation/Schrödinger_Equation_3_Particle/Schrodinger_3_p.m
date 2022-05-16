clc
clear all

N = 100;

m       = 1;
hbar    = 1;
alpha   = -4;
bound   = 5;
x       = linspace(-bound, bound,N);
x       = [x x];
y       = x;
dx      = x(2) - x(1);

mtx = sparse(N^2, N^2);
a = 1;
b = N;
for i = 1:2
    for j = a:b
        mtx(j,j) = 0.25 * (x(j)^2 - alpha)^2;
    end
    a = a + N;
    b = b + N;
end

for i = 2:N^2
    mtx(i - 1, i) = -1;
    mtx(i, i - 1) = -1;
end

eta = 20;
for i = 1:N
    Interaction(i
end
a = eig(mtx);
plot(a)