clc
clear all
format long

%% Structural Initialization:
FS  = 1;        % Figure switch
FE  = 1;        % ErrorFinding Figure switch
PN  = 3;        % Particle Number

NoA = 4001;
Eta = 20;
%% Loadig EQPOS data
addpath('D:\BME PhD\.Wigner Crystal Collective Tunneling\CollectiveTunneling\1 - Equilibrium Positions\Data')
if PN == 1 % Struct
    AlphaValues = linspace(-3, 3, NoA);
    for i = 1:length(AlphaValues)
        if AlphaValues(i) < 0
            EqPos(i) = 0;
        else
            EqPos(i) = sqrt(AlphaValues(i));
        end
    end
    EqPos = EqPos';
elseif PN == 3 % Struct
    Data        = load('ThreeParticleEquilibriumPositions.mat');
    AlphaValues = Data.DataStruct.Alpha;
    EqPos       = Data.DataStruct.EquilibriumPositions;
elseif PN == 5 % Struct
    Data        = load('FiveParticleEquilibriumPositions.mat');
    AlphaValues = Data.DataStruct.Alpha;
    EqPos       = Data.DataStruct.EquilibriumPositions;
elseif PN == 7 % Matrix
    Data        = load('SevenParticleEquilibriumPositions.mat');
    AlphaValues = -Data.eqpos(end, :);
    EqPos       = Data.eqpos(1:end-1, :)';
end

if FS == 1
    figure(1)
    clf(figure(1))
    hold on
    plot(AlphaValues, EqPos)
    yline(0)
    hold off
end

%% Calculating Freq. and Modes vecs

for alphaInd = 1:length(AlphaValues)
    V           = [];
    U_diag      = zeros(PN, 1);
    U_offdiag   = zeros(PN, PN);

    for particleIndI = 1:PN
        V(particleIndI) = (3 * EqPos(alphaInd, particleIndI)^2 - AlphaValues(alphaInd));

        for particleIndJ = 1:PN
            if particleIndI == particleIndJ

            else
                U_diag(particleIndI)                    = U_diag(particleIndI)  + (Eta * (2/abs((EqPos(alphaInd, particleIndI) - EqPos(alphaInd, particleIndJ))^3)));
                U_offdiag(particleIndI, particleIndJ)   =                       - (Eta * (2/abs((EqPos(alphaInd, particleIndI) - EqPos(alphaInd, particleIndJ))^3)));
            end
        end
    end
    
    for i = 1:PN
        for j = 1:PN
            if i == j
                SpringMatrix(i, j) = V(i) + U_diag(i);
            else
                SpringMatrix(i, j) = U_offdiag(i, j);
            end            
        end
    end

    % Egzakt diagonalization:
    [ModVec, EigVal]            = eig(SpringMatrix);
    Frequencies(alphaInd, :)    = diag(EigVal);
    EigModes(alphaInd, :)    = ModVec(:, 1);
end
%%
if FS == 1
    figure(19)
    clf(figure(19))
    hold on
    for n = 1:PN
        CurveName = [num2str(n) ' \omega'];
        plot(AlphaValues, sqrt(Frequencies(:, n)), '.-', 'DisplayName', CurveName)
    end
    minimum = min(min(Frequencies));
    [x, y] = find(Frequencies == minimum);
    xline(AlphaValues(x), 'DisplayName', ['\alpha_c \sim ' num2str(AlphaValues(x))])
    %grid on
    %legend
    xlabel('\alpha', 'FontSize', 40)
    ylabel('\omega_i', 'FontSize', 40)
    set(gca,'fontsize', 30)
    xlabel('\alpha', 'FontSize', 40)
    ylabel('\omega_i', 'FontSize', 40)
    hold off

end
%%
if FS == 1
    figure(3)
    clf(figure(3))
    hold on
    for n = 1:PN
        plot(AlphaValues, abs(EigModes(:, n)), '.-', 'DisplayName', 'ModVec comp.')
    end
    grid on
    legend
    xlabel('\alpha', 'FontSize', 22)
    ylabel('\omega_i', 'FontSize', 22)
    hold off

end


FileName = ['Omega_sq_Freqs_N_' num2str(PN)];
save(FileName, 'EigModes')







