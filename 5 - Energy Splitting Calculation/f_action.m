function [propagator, action, splitt1, splitt2, splitt3, splitt4, splitt5] = f_action(eta, alpha, trajectory, r, z_time, EigValu, propagator, S, VS, analytical)
%first the action shift:
shift = f_initshift(eta, alpha, trajectory);

disp(shift)
%then the actual action of the trajectory:
action = f_actioncalc(trajectory, r, alpha, eta, 200, z_time, (z_time(2) - z_time(1)), shift);
z_time = linspace(-1 + 10^-15 , 1 - 10^-15, length(S));
z_time = linspace(-1, 1, length(S));

disp([' action = ', num2str(action)])

omega = sqrt(EigValu(1,1));
[gof, fc]   = f_fitting_VS(S, VS);

VS2         = fc.b * (S.^2 - min(S)^2).^2;
%VS2         = fc.b * (S.^2 - 0.5 * (fc.c - min(S)^2)).^2;
VS = VS - min(VS);
if analytical == 1
    VS2 = VS2 - min(VS2);
end
figure(10)
clf(figure(10))
hold on
plot(S, VS - min(VS), '.-')
plot(S, VS2 - min(VS2), 'o-')

% plot(S, VS - VS2)
hold off

%omegaS_sq       = 4 * fc.b * (2 * fc.c);
omegaS_sq       = 4 * fc.b * (2 * min(S)^2);
%omega           = sqrt(omegaS_sq);
%Easy Milnikov:
splitt1 = 4 * omega * max(S) * sqrt(omega/pi) * exp(-action);%sqrt(2)^2 * (omega) * sqrt(abs(alpha)) * sqrt((omega))/sqrt(pi) * exp(-(action));
%Landau approx.:
splitt2 = omega/pi * exp(-action);
%Instanton pref.:
splitt3 = 2 * sqrt((6 * action)/pi) * omega * exp(-action);

%Proper Milnikov:----------------------------------------------------------
% p1      = sqrt(2 * 0.25 .* (trajectory(1, :).^2 + alpha).^2);
% p2      = sqrt(2 * 0.25 .* (trajectory(2, :).^2 + alpha).^2);
% p3      = sqrt(2 * 0.25 .* (trajectory(3, :).^2 + alpha).^2);
% pA      = p1 + p2 + p3;
% 
% for i = 2:length(trajectory)/2
%     dp1(i) = (p1(i) - p1(i - 1))/abs(trajectory(1, i) - trajectory(1, i - 1));
%     dp2(i) = (p2(i) - p2(i - 1))/abs(trajectory(2, i) - trajectory(2, i - 1));
%     dp3(i) = (p3(i) - p3(i - 1))/abs(trajectory(3, i) - trajectory(3, i - 1));
%     
%     dpA(i) = (pA(i) - pA(i - 1))/(S(i) - S(i - 1));
% 
%     PPP1(i) = (sqrt(dp1(end)^2 + dp2(end)^2 + dp3(end)^2));
%     PPP2(i) = abs(dpA(end));
% 
% end
% func = (omega - PPP1)./(1 - z_time(2:end/2+1).^2);
% func2 = (omega - PPP2)./(1 - z_time(2:end/2+1).^2);
% func(2) = 0;
% % func(1) = 0;
% figure(2)
% clf(figure(2))
% hold on
% plot(func)
% hold off
% 
% int = 0;
% int2= 0;
% for i = 3:100
%     int = int + (func(i) - func(i-1)) * (z_time(i) - z_time(i-1));
%     int2= int2+ (func2(i) - func2(i-1)) * (z_time(i) - z_time(i-1));
% end
% int
% int2
% splitt4 = exp((int2)) * exp(-action) * sqrt(omega/pi) * 2 * omega * abs(alpha); %abs(trajectory(2, 1));

P_s = sqrt(2 * VS);

dP_s(1) = (P_s(3) - P_s(1))/(S(3) - S(1));
for i = 2:length(P_s)-1
    dP_s(i) = (P_s(i + 1) - P_s(i - 1))/(S(i + 1) - S(i - 1));
end
dP_s(end+1) = (P_s(end) - P_s(end - 2))/(S(end) - S(end - 2));

% if analytical == 1
%     for i = 1:length(P_s)
%         dP_s(i) = sqrt(2) *  2 * fc.b * S(i) * (S(i)^2 - fc.c)/sqrt(fc.a + fc.b * (S(i)^2 - fc.c)^2);
%     end
% end

figure(100)
clf(figure(100))
hold on
%dP_s = smooth(dP_s);
plot(S, dP_s, '.')
hold off

%Integral:
int     = 0;
int2    = 0;
func    = (r./(1 - z_time.^2)) .* (omega - dP_s);
func2   = (r./(1 - z_time.^2)) .* (dP_s(1) - dP_s);
% func(1) = 0;
% func(end) = 0;
figure(100)
clf(figure(100))
hold on
%dP_s = smooth(dP_s);
plot(z_time, func, '.')
plot(z_time, func2, 'o')
hold off
xlim([-1 0])

for i = 2:length(P_s)/2 - 1
    dz = (z_time(i + 1) - z_time(i));
    %dS = S(i + 1) - S(i);
    int = int +  dz * 0.5 * (func(i + 1) + func(i));
    int2 = int2 + dz * 0.5 * (func2(i + 1) + func2(i));
end
disp([' Integral: ',num2str(int)])
disp([' Integral2: ',num2str(int2)])
splitt4 = sqrt(4 * omega) / pi * max(P_s(1:end/2-1)) *  exp(int) * exp(-action);
splitt5 = sqrt(4 * omega) / pi * max(P_s) *  exp(int2) * exp(-action);
disp(['1D part= ', num2str(splitt1)])

propagator = propagator;
end

