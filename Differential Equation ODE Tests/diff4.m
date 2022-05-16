function [timeval,yval] = diff4(a, y0, z0, time)
for i = 1:length(a)-1
    if i == 1
        init(1) = y0;
        init(2) = z0;
    else
        init(1) = y(end,1);
        init(2) = y(end,2);
    end
    param = a(i);
    tspan = [time(i) time(i+1)];
    
    [t,y] = ode45(@(t,y)[param^2 - y(1)^2; y(1) - y(2)], tspan, [init(1) init(2)]);
    if i == 1
        timeval     = t';
        yval        = y';
    else
        timeval     = [timeval t'];
        yval        = [yval y'];
    end
end
end

