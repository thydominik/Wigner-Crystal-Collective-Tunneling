function [time, yval] = diff1(a, ts, tf, y0)
    [time, yval] = ode45(@(t,y) [a^2 - y^2], [ts tf], y0);
end 

