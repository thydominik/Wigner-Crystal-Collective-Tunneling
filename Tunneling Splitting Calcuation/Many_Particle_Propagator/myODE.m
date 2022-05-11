function [dydt] = myODE(a,b,c,d,e,f,g,h,i, omega)
dydt        = zeros(9,1);

dydt(1) = omega(1,1) - xi(1)^2 - xi(2)*xi(4) - xi(3)*xi(7);
dydt(2) = omega(1,2) - xi(1)*xi(2) - xi(2)*xi(5) - xi(3)*xi(8);
dydt(3) = omega(1,3) - xi(1)*xi(3) - xi(2)*xi(6) - xi(3)*xi(9);
dydt(4) = omega(2,1) - xi(1)*xi(4) - xi(4)*xi(5) - xi(6)*xi(7);
dydt(5) = omega(2,2) - xi(2)*xi(4) - xi(5)^2 - xi(6)*xi(8);
dydt(6) = omega(2,3) - xi(3)*xi(4) - xi(6)*xi(5) - xi(6)*xi(9);
dydt(7) = omega(3,1) - xi(1)*xi(7) - xi(4)*xi(8) - xi(7)*xi(9);
dydt(8) = omega(3,2) - xi(2)*xi(7) - xi(8)*xi(5) - xi(8)*xi(9);
dydt(9) = omega(3,3) - xi(3)*xi(7) - xi(6)*xi(8) - xi(9)^2;

% dydt(1) = omega(1,1) - a^2 - b*d - c*g;
% dydt(2) = omega(1,2) - a*b - b*e - c*h;
% dydt(3) = omega(1,3) - a*c - b*f - c*i;
% dydt(4) = omega(2,1) - a*d - d*e - f*g;
% dydt(5) = omega(2,2) - b*d - e^2 - f*h;
% dydt(6) = omega(2,3) - c*d - f*e - f*i;
% dydt(7) = omega(3,1) - a*g - d*h - f*i;
% dydt(8) = omega(3,2) - b*g - h*e - h*i;
% dydt(9) = omega(3,3) - c*g - f*h - i^2;
end


