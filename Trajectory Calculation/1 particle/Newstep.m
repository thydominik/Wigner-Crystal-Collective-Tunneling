function Pos = newstep(Pos, NoP, sigma, alpha, z)
    R = randi((NoP/2)-1) + 1;
    D = (normrnd(Pos(R),sigma));
    temp = Pos(R);
    Pos(R) = D;
    while Pos(R) > 0 || Pos(R) < -sqrt(alpha) || Pos(R) > (z(R) * abs(Pos(1)))
       D = (normrnd(temp,sigma));
       Pos(R) = D;
    end
    Pos(NoP-(R-1)) = -Pos(R);

end