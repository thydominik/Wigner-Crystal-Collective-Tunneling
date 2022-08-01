function [khi, shift] = initpos(a, NoP)
    khi = linspace(-sqrt(a), sqrt(a), NoP);
    shift = 0;
end