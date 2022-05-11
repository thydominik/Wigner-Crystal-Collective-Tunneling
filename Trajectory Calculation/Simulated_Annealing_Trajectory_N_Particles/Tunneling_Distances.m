clc
clear all

%Loading the equlibrium positions:
eqpos   = load('EqPos_eta20_alpha_5_20_5p.mat');
eq_pos  = eqpos.eqpos;

%Getting the number of particles involved:
N = length(eq_pos(:, 1)) - 1;

%Equilibrium positions:
for State = 1:length(eq_pos)
    for q = 1:N
        p_in(q)     = eq_pos(q, State);
        p_out(q)    = -eq_pos(N - q + 1, State);
        dist(State, q) = p_in(q) - p_out(q);
    end
end

figure(1)
clf(figure(1))
hold on
title('Tunneling distances')
plot(-eq_pos(end, :), abs(dist(1:end, 1)), '.', 'DisplayName', '1st')
plot(-eq_pos(end, :), abs(dist(1:end, 2)), '.', 'DisplayName', '2nd')
plot(-eq_pos(end, :), abs(dist(1:end, 3)), '.', 'DisplayName', '3rd')
legend
grid on
hold off
