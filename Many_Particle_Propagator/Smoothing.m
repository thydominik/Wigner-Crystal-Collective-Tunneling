function [smooth_trajectory, N_division, coefficients] = Smoothing(trajectory, ytime, lim, N, eq_pos, Tau0)   

y_data1 = trajectory(1,:);
y_data2 = trajectory(2,:);
x0 = eq_pos;
[cf1] = createfit1(ytime, y_data1, x0(1), -x0(3));
[cf2] = createfit1(ytime, y_data2, x0(2), -x0(2));
 
x = linspace(-lim, lim, N);
 
smooth_trajectory(1,:) = cf1(1) + cf1(2) * tanh(x*cf1(3) + cf1(4));
smooth_trajectory(2,:) = cf2(1) + cf2(2) * tanh(x*cf2(3) + cf2(4));
smooth_trajectory(3,:) = -cf1(1) + cf1(2) * tanh(x*cf1(3) - (cf1(4)));
N_division = N;

%cf1(4) = cf1(4) - cf2(4);
%cf2(4) = 0;

cf3 = [-cf1(1) cf1(2) cf1(3) -(cf1(4))];
coefficients = [cf1' cf2' cf3'];
disp(num2str(coefficients))
end

