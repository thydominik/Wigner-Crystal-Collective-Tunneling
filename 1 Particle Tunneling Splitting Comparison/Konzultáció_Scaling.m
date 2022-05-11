clc
clear all
scale = 1;
shift = -156.13;
V = [110 330 610 850];
dV = [0 230 330 280]/scale;
dac = [0 4.5 3.3 2.7];
ac = [0 -4.5 -7.8 -10.5];
aq = [-3 -7.2 -10.8 -13.5];
N = [1 3 5 7];
plot(N, dV)

figure(1)
clf(figure(1))
hold on
%plot(V, ac, '.')
scale = (V(3)-shift)/(-aq(3))
plot((V - shift)/scale, -aq, 'o-')
ylim([0 15])
xlim([0 15])

hold off

figure(2)
clf(figure(2))
hold on
plot(N,(V - shift)/scale, 'o-', 'DisplayName', '(V-V0)/scale')
plot(N, -aq,' .-', 'DisplayName', '\alpha_q')
xlim([0 8])
ylim([1 15])
legend
hold off