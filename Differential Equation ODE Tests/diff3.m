function [tval,yval,zval] = diff3(a,ts,tf,y0,z0)
    [tval,yval,zval] = ode45(@(t,y)[a^2 - y(1)^2; y(1) - y(2)],[ts tf],[y0 z0]);
end

