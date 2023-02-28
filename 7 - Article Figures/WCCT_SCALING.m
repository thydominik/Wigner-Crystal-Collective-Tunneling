clc
clear all

V_kink      = [110 322.5 640 935];
x_n         = [23.4 39.1 55.5 65.6];
alpha_shift = 2;
yQ          = [90 13 7.5, 7.5];


% V* = alpha_shift + alpha_ED * x_n 

alpha_kink = (V_kink)./ x_n - alpha_shift