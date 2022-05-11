function [Positions, Discarded, Energy, R_param] = TrajsMC(eq_pos, eta, Alpha, TrajsIter, Z, dz, N)
%TRAJSMC

    %OUTPUT ------------------------------------------------
    %
    %OUTPUT ------------------------------------------------
    
    %INPUT ARGUMENTS ---------------------------------------
    %
    %INPUT ARGUMENTS ---------------------------------------
    
    Temp_init   = 50;
    Temp        = Temp_init * exp(-linspace(0, 30, TrajsIter));
    Sigma       = 0.1 * sqrt(Temp);
    
    %this r value needs to be changed with alpha, but I have no idea how
    %so, this will be a guess:
    R_param = linspace(1.3, 0.8, length(Alpha));
    
    Energy      = [];
    Discarded   = [];
    
    Positions = zeros(length(Alpha), 3, N);
    for i = 1:length(Alpha)
        r = R_param(i);
        
        p1_in =  eq_pos(1, i);        %initial and final positions of the particles
        p1_fi = -eq_pos(3, i);
        p2_in =  eq_pos(2, i);
        p2_fi = -eq_pos(2, i);
        p3_in =  eq_pos(3, i);
        p3_fi = -eq_pos(1, i);
        
        [position, shift]   = f_initpos(N, p1_in, p1_fi, p2_in, p2_fi, p3_in, p3_fi, eta, Alpha(i));
        E_0                 = f_actioncalc(position, r, Alpha(i), eta, N, Z, dz, shift);
        disp("E_0 = " + num2str(E_0, 15))
        disp("Shift = " + num2str(shift))
        
        for j = 1:TrajsIter

            sig         = Sigma(j);
            pos_new     = f_newstep(position, N, sig, p1_in, p1_fi, p2_in, p2_fi, p3_in, p3_fi);
            E_new       = f_actioncalc(pos_new, r, Alpha(i), eta, N, Z, dz, shift);
            E_diff      = E_0 - E_new;

            if E_diff > 0
                position    = pos_new;
                E_0         = E_new;
            elseif rand() <= exp(E_diff/Temp(j))
                position    = pos_new;
                E_0         = E_new;
            else
                %Discarded(j, i) = j;
            end
            %Energy(j, i) = E_0;

            if rem(j,20000) == 0
                figure(1)
                clf(figure(1))
                hold on
                plot(Z,position(1,:))
                plot(Z,position(2,:))
                plot(Z,position(3,:))
                hold off
                disp("iter= " + num2str(j) + "   "+ "E_0= " + num2str(E_0, 5))
            end   
        end
        Energy(i) = E_0;
        Positions(i, :, :) = position;
        disp(['done with alpha = '  num2str(Alpha(i)) ' ' num2str(Energy(end, i))])
        
    end
    
end

