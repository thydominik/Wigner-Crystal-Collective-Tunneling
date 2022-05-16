function [shift] = f_initshift(rs, a, pos)  
    in1 = pos(1,1);
    in2 = pos(2,1);
    in3 = pos(3,1);
    
    F1 = 0.25 *(in1^2 + a)^2;
    F2 = 0.25 *(in2^2 + a)^2;
    F3 = 0.25 *(in3^2 + a)^2;
    
    Q = rs * ((1/abs(in2 - in1)) + 1/abs(in3 - in1) + 1/abs(in3 - in2 )) ;
    shift =  Q + (F1 + F2 + F3);

end

