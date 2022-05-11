function [pos, shift] = f_initpos(N, in1, fin1, in2, fin2, in3, fin3, rs, a)
    pos = zeros(3,N);
    
    for i = 1:N/2 
        pos(1,i) = in1;%abs(in1 - fin1)*rand() + in1;
        pos(2,i) = in2;%abs(in2 - fin2)*rand() + in2;
        pos(3,i) = in3;%abs(in3 - fin3)*rand() + in3;
    end
    for i = N/2+1:(N)
        pos(1,i) = fin1;%abs(in1 - fin1)*rand() + in1;
        pos(2,i) = fin2;%abs(in2 - fin2)*rand() + in2;
        pos(3,i) = fin3;%abs(in3 - fin3)*rand() + in3;
    end
    
    pos(1,1) = in1;
    pos(1,N) = fin1;
    pos(2,1) = in2;
    pos(2,N) = fin2;
    pos(3,1) = in3;
    pos(3,N) = fin3;
    
    F1 = 0.25 *(pos(1,1)^2 + a)^2;
    F2 = 0.25 *(pos(2,1)^2 + a)^2;
    F3 = 0.25 *(pos(3,1)^2 + a)^2;
    
    Q = rs * ((1/abs(pos(2,1) - pos(1,1))) + 1/abs(pos(3,1)-pos(1,1)) + 1/abs(pos(3,1)-pos(2,1))) ;
    shift =  Q + (F1 + F2 + F3);
end

