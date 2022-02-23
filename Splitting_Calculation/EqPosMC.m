function [eq_pos, Discarded, Final_Energy] = EqPosMC(PartNum, Alpha, eta, EqPosIter )
%EQPOSMC: Does the Monte Carlo 

    %OUTPUT ------------------------------------------------
    %
    %OUTPUT ------------------------------------------------
    
    %INPUT ARGUMENTS ---------------------------------------
    %
    %INPUT ARGUMENTS ---------------------------------------
    
    Temp_init   = 600;                                              %starting temperature
    Temp        = Temp_init * exp(-(linspace(0,80,EqPosIter))/2);   %exponential temperature decrease  
    Sigma       = 0.1 * sqrt(Temp);
    
    %for loop in the alpha vector
    for i = 1:length(Alpha)
        Positions   = initpos(PartNum);
        Energy_0    = energy(Positions, PartNum, Alpha(i), eta);
        
        Discarded   = [];
        
        %simulated annealing steps
        for j = 1:EqPosIter
            New_Position    = new_pos(Positions, Temp(j), Alpha(i), eta, Sigma(j));
            Energy          = energy(New_Position, PartNum, Alpha(i), eta);
            E_diff          = Energy_0 - Energy;

            if E_diff > 0
                Positions = New_Position;
                Energy_0 = Energy;
            elseif rand() < exp(E_diff/Temp(j))
                Positions = New_Position;
                Energy_0 = Energy;
            else
                %Discarded(j, i) = i;
            end
     
        end
        Final_Energy(i)   = Energy_0;
        eq_pos(:,i) = Positions;
        
        disp(['done with alpha = ' num2str(Alpha(i))])
    end

    %ordering the positions to one side so that the left side has 2
    %particles and the right has only 1 always:
    eq_pos = organizer(eq_pos);
    
end

