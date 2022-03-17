clc
clear all

alpha   = -11.9;
eta     = 20;

trajectory  = load('Traj_5p_70.mat');
trajectory = trajectory.Position;
LS(1)       = 0;
for i = 2:length(trajectory)
    r(1) = trajectory(1, i) - trajectory(1, i - 1);
    r(2) = trajectory(2, i) - trajectory(2, i - 1);
    r(3) = trajectory(3, i) - trajectory(3, i - 1);
    r(4) = trajectory(4, i) - trajectory(4, i - 1);
    r(5) = trajectory(5, i) - trajectory(5, i - 1);
    S(i)    = sqrt(r(1)^2 + r(2)^2 + r(3)^2 + r(4)^2 + r(5)^2);
    LS(i)   = LS(i - 1) + S(i);
end

for i = 1:length(trajectory)
    V = 0;
    Q = 0;
    for j = 1:5
        V = V + 0.25 * (trajectory(j, i).^2 + alpha).^2;
    end
    for j = 1:5
        for k = (j+1):5
            Q = Q + eta/abs(trajectory(j, i) - trajectory(k, i));
        end
    end
    
    VS(i) = V + Q;

end


figure(1)
clf(figure(1))
VS_S = VS - min(VS);
VS_r = VS_S(10:end-10);
S_S = LS - max(LS)/2;
S_r = S_S(10:end-10);
x = linspace(-1, 1, 100);
weight = (1 - exp(-0.5 * x.^2))/0.02;
w_r = weight(10:end-10);
hold on
x = linspace(S_S(1), S_S(end), 100);
x = S_S;
a = max(VS_S);
b = -3.722;
c = 0.4615;
func = a + b * x.^2 + c * x.^4;

plot(S_S(10:end-10), (VS_S(10:end-10)), 'r-', 'LineWidth', 2, 'DisplayName', 'Effective potential (5 part.)')
plot(x(10:end-10),  func(10:end-10), 'ko', 'LineWidth', 0.9, 'DisplayName', '4th order polynomial fit')

xlabel('S', 'FontSize', 12)
ylabel('V(S)', 'FontSize', 12)
xlim([-2.2 2.2])
ylim([0 8])
legend
hold off

%%

figure(2)
clf(figure(2))
hold on
z = linspace(-1, 1, 100);

plot(z, smooth(trajectory(1, :)), '-', 'LineWidth', 2)
plot(z, smooth(trajectory(2, :)), '-', 'LineWidth', 2)
plot(z, smooth(trajectory(3, :)), '-', 'LineWidth', 2)
plot(z, smooth(trajectory(4, :)), '-', 'LineWidth', 2)
plot(z, smooth(trajectory(5, :)), '-', 'LineWidth', 2)
%legend
xline(0)
yline(0)
xlabel('z', 'FontSize', 12)
ylabel('\chi_i', 'FontSize', 18)
hold off
%%
figure(3)
clf(figure(3))
hold on
plot(z, smooth(smooth(LS)), 'k-', 'LineWidth', 2)
xlabel('z', 'FontSize', 12)
ylabel('S(z)', 'FontSize', 12)
hold off