function [timeval, yval] = diff2(time, a, y0)
    for i = 1:length(a)-1
        if i == 1
            init = y0;
        else
            init = y(end);
        end
        param = a(i);
        tspan = [time(i) time(i+1)];

        [t,y] = ode45(@(t,y)[param^2 - y^2], tspan, init);
        if i == 1
            timeval = t';
            yval = y';
        else
            timeval = [timeval t'];
            yval = [yval y'];
    end
    
    end
end

