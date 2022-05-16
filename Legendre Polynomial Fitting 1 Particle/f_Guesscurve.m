function [C_param, q, dQ] = f_Guesscurve(alpha, N, div, r, z)
    C = zeros(N, 1);
    q = zeros(N, div);
    %let's construct the least complicated initial guess curve 
    %C(1) = sqrt(alpha); %the second Legendre polynomial is 'L_1(x) = x'
    C_param = C';
    
    q = sqrt(alpha) * tanh(atanh(z)*sqrt((alpha + 0.3)/2)*r);
    Q = q;
    dQ = zeros(1,N);
    for i = 1:div
        prefactor = sech(sqrt(alpha/2) * r *atanh(z(i)))^2;
        if prefactor == 0
            dQ(i) = 0;
        else
            dQ(i) = prefactor*(r/(1-z(i)^2))*alpha /sqrt(2) * r;
        end
    end
    
end

