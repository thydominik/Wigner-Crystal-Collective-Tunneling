function [C_param, Q] = guesscurve(alpha, N, div, Legendre)
    C = zeros(N, 1);
    q = zeros(N, div);
    %let's construct the least complicated initial guess curve 
    C(1) = sqrt(alpha); %the second Legendre polynomial is 'L_1(x) = x'
    C_param = C';
    for i = 1:N
        q(i,:) = q(i,:) + C(i) * Legendre(i,:);
    end
    Q = sum(q);
end

