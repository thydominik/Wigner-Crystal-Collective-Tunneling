function propagator = f_action(eta, alpha, trajectory, r, z_time, EigValu, propagator)
%first the action shift:
shift = f_initshift(eta, alpha, trajectory);

%then the actual action of the trajectory:
action = f_actioncalc(trajectory, r, alpha, eta, 3, z_time, (z_time(2) - z_time(1)), shift);

action = 4.653598741;

disp([' action = ', num2str(action)])

omega = sqrt(EigValu(1,1));
splitt = sqrt(2)^2 * omega * sqrt(abs(alpha)) * sqrt(omega/pi) * exp(-action);
disp(['1D part= ', num2str(splitt)])

propagator = propagator * splitt;
end

