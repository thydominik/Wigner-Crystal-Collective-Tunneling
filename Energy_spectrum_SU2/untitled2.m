clc
clear all

M = [1 2 5; 0 2 7; 0 0 5]

[V, E] = eig(M)

D = E
P = [V(:, 1) V(:, 2) V(:, 3)]
Pi = inv(P)
P*D*Pi 